#include <boost/foreach.hpp>
#include <stdint.h>
#include <glob.h>
#include <pnetcdf.h>
#include <string>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>
#include "app.h"
#include "bil/bil.h"
#include "tinyxml/tinyxml.h"
#include "partitioner/partitioner.h"
#include "common/utils.h"
#include "cons_kdtree.h" 

// void forbp(){
//   fprintf(stderr, "In ForBP\n");
// }

#define WRITE 1

typedef diy::DiscreteBounds DBounds; 
typedef diy::RegularGridLink RGLink;
typedef diy::RegularDecomposer<DBounds> Decomposer;

struct AddBlock { // functor
  AddBlock(diy::Master &master_) :
    master(master_)
  {}
  
  void operator()(int gid, 
      const DBounds& core, 
      const DBounds& bounds, 
      const DBounds& domain, 
      const RGLink& link) const {


    Block *b = new Block();
    RGLink *l = new RGLink(link);
    diy::Master &m = const_cast<diy::Master&>(master);

    int lid = m.add(gid, b, l);
    b->gid = gid;
    for (int i=0; i<4; i++) {  // redefined in the first k-d tree

      b->lb[i] = l->core().min[i];
      b->ub[i] = l->core().max[i];
      b->glb[i] = l->bounds().min[i];
      b->gub[i] = l->bounds().max[i];
      b->lload[i] = l->bounds().min[i];
      b->uload[i] = l->bounds().max[i];
    }
    b->data_bounds = domain;
  }

  diy::Master& master;
};

CPTApp::CPTApp() 
  : _comm_world(MPI_COMM_WORLD),
    _time_io(0),
    _nt_loops(0)
{
  for (int i=0; i<4; i++) {
    _domain_size[i] = 1;
    _ghost_size[i] = 1;

    _block_size[i] = 1;
  }
}

CPTApp::~CPTApp()
{
  delete _decomposer;
  delete _assigner;
  delete _master;
  delete _communicator;
}

void CPTApp::init(int argc, char **argv)
{
  _timestamp_start = MPI_Wtime();
  init_mpi(); 

   // if (_comm_world_rank == 1)
   //  {
        

   //      volatile  int dwait=0;
   //      fprintf(stderr , "pid %ld  waiting  for  debugger\n"
   //          , (long)getpid ());
   //          while(dwait==0) { /*  change  ’i’ in the  debugger  */ }

   //  }
   //  MPI_Barrier(_comm_world);


  if (_comm_world_rank == 0)
    parse_arguments(argc, argv);
  bcast_config();

  int num_dims = space_only() ? 3 : _num_dims;
  if (_appconf.block_mem_limit() < 0) _constrained  = false;
  else _constrained = true;

  // determine unified ghost size by user-defined partitioner
  if (_appconf.elastic_ghost()) {
    probe_elastic_ghosts(_block_size, _ghost_size);
    //probe_elastic_ghost_block_size(1, _block_size, _ghost_block_size); // grid size
  } else {
    probe_bs(_block_size);
    for (int i=0; i<num_dims; i++) 
      _ghost_size[i] = _appconf.ghost_size();
  }

  _communicator = new diy::mpi::communicator(_comm_world);
  _master = new diy::Master(*_communicator, _num_threads);
  _assigner = new diy::RoundRobinAssigner(_communicator->size(), _num_blocks);

  for (int i = 0; i < num_dims; i ++) { 
    _domain.min[i] = 0;
    _domain.max[i] = _domain_size[i] - 1;
  }

  diy::RegularDecomposer<DBounds>::BoolVector share_face;
  diy::RegularDecomposer<DBounds>::BoolVector wrap;

  diy::RegularDecomposer<DBounds>::CoordinateVector ghosts;
  for (int i = 0; i < num_dims; i ++)
    ghosts.push_back(_ghost_size[i]); // ghosts 
  _decomposer = new diy::RegularDecomposer<DBounds>(
      num_dims, 
      _domain, 
      //*_assigner,
      _num_blocks, 
      share_face,
      wrap, 
      ghosts); 
  _divisions = _decomposer->get_divisions();

  if (comm_world_rank() == 0) 
  {
    fprintf(stderr, "================ INIT ==============\n");
    if (is_kd_tree()) 
      fprintf(stderr, "%s, with k-d tree decomposition\n", _dataconf.name().c_str());
    else
      fprintf(stderr, "%s, without k-d tree decomposition\n", _dataconf.name().c_str());
    fprintf(stderr, "rk=%d, np=%d, nt=%d,\n", comm_world_rank(), comm_world_size(), num_threads());
    fprintf(stderr, "block_memory_limit=%d, space_only=%d\n", _appconf.block_mem_limit(), space_only() ? 1 : 0);
    fprintf(stderr, "divisions=");
    for(int i=0; i<_divisions.size(); i++)
      fprintf(stderr, "%d,", _divisions[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "nb=%d, ds={%d,%d,%d,%d}, bsz={%d,%d,%d,%d}, gs={%d,%d,%d,%d},\n", _num_blocks, 
      _domain_size[0], _domain_size[1], _domain_size[2], _domain_size[3], 
      _block_size[0], _block_size[1], _block_size[2], _block_size[3], 
      _ghost_size[0], _ghost_size[1], _ghost_size[2], _ghost_size[3]);
    fprintf(stderr, "====================================\n");
  }

  AddBlock addblock(*_master);
  //if (_appconf.elastic_ghost()) {
  //  _decomposer->set_elastic_ghosts(true);
    //_decomposer->set_ghosts(_ghosts);
  //} else 
  //  _decomposer->set_elastic_ghosts(false); 
  _decomposer->decompose(_communicator->rank(), *_assigner, addblock);
  _assigner->local_gids(comm_world_rank(), _gids); 
  //if (!is_baseline() && is_kd_tree())  // k-d tree to reassign bounds to blocks
  if (is_kd_tree() && _constrained)
    pt_cons_kdtree_exchange(*_master, *_assigner, _divisions, num_dims, space_only(), _block_size, _ghost_size, true, true, false);


   
  std::vector<std::string> local_entities(gids().size());
  int n = 0;
  MPI_Barrier(_comm_world);
  double iostart = MPI_Wtime();
  _timestamps.push_back(iostart); // start
  BIL_Init(_comm_world); // BIL is used with comm_world
  // fprintf(stderr, "!!!ASDF10!!!!\n");
  
  BOOST_FOREACH(int gid, gids()) { // may be multiple blocks per proc. if _nb_per_proc > 1
    Block &b = block(gid);

    float lpoint[num_dims];
    for (int i = 0; i < num_dims; i ++) lpoint[i] = b.lb[i];
    int bound_gid = pt2gid(lpoint);
    //bound_gids.push_back(pt2gid(lpoint));

    PBBoundGid bg;
    bg.set_gid(gid);
    bg.set_bound_gid(bound_gid);
    bg.SerializeToString(&local_entities[n++]);

    if (space_only()) {
      b.lb[3] = 0;
      b.ub[3] = _domain_size[3] - 1;
      b.glb[3] = 0;
      b.gub[3] = _domain_size[3] - 1;
      b.lload[3] = 0;
      b.uload[3] = _domain_size[3] - 1;
    }
//b.print_info();
    if (_dataconf.name() == "nek5000" && 0) {
      b.add_block_nek(_dataconf);
      //_time_io += iotime;
      //fprintf(stderr, "rank %d read finished.\n", comm_world_rank());
    } else {
//size_t mem = 3*4*(b.uload[3] - b.lload[3] + 1) * (b.uload[2] - b.lload[2] + 1) * (b.uload[1] - b.lload[1] + 1) * (b.uload[0] - b.lload[0] + 1);
//fprintf(stderr, "block mem = %f\n", mem/1024.f/1024.f);
      if (_num_dims == 4) b.bil_add_block_4D(_dataconf);
      else if (_num_dims == 3) b.bil_add_block_3D(_dataconf);
      else if (_num_dims == 2) b.bil_add_block_2D(_dataconf);
      else assert(false);
    }
  }



  BIL_Read();
  BIL_Finalize();

  MPI_Barrier(_comm_world);
  double ioend = MPI_Wtime();
  _time_io += ioend - iostart;

#if 0
  // gid = proc. TODO if multiple blocks per process
  _bound_gids.resize(comm_world_size()); 
  for (int i = 0; i < comm_world_size(); i ++) _bound_gids[i] = i;
  //if (is_kd_tree()) // TODO improve efficiency
  {
    std::vector<int> each_num(comm_world_size());
    MPI_Allgather(&bound_gid, 1, MPI_INT, each_num.data(), 1, MPI_INT, comm_world());
    for (int i = 0; i < comm_world_size(); i ++) {
      for (int j = 0; j < comm_world_size(); j ++) {
        if (each_num[j] == i) {
          _bound_gids[i] = j; 
          break;
        }
      }
    }
  }
#endif 

  _bound_gids.resize(_num_blocks);
  int local_num = gids().size(), total_num;
  std::vector<std::string> all_entities;
  gather_protobuf(local_num, local_entities, total_num, all_entities, -1, comm_world());
  for (int i = 0; i < total_num; i ++) {
    PBBoundGid bg;
    bg.ParseFromString(all_entities[i]);
    _bound_gids[bg.bound_gid()] = bg.gid();
  }

  /*
  if (comm_world_rank() == 0) {
    fprintf(stderr, "bound_gids: ");
    for (int i=0; i<_bound_gids.size(); i++)
      fprintf(stderr, "%d, ", _bound_gids[i]);
    fprintf(stderr, "\n");
  }
  */
  MPI_Barrier(_comm_world);
  _timestamp_init = MPI_Wtime();
  _timestamps.push_back(_timestamp_init); // io
  _timecategories.push_back(0);


}

void CPTApp::probe_elastic_ghost_block_size(
    const int factor,  
    int *bsz, int *gsz
  ) const // old refined by Jiang
{
  const int64_t mem_limit = (int64_t)_appconf.block_mem_limit() * 1024 * 1024 / (_num_blocks / comm_world_size());

  int num_dims = space_only() ? 3 : _num_dims;

  RegularPartitioner partitioner(
      _comm_world, num_dims, _dataconf.runs_size(), _dataconf.variables_size());
  partitioner.set_domain_size(_domain_size);
  partitioner.set_num_blocks(_num_blocks);
  int block_start[4], block_size[4];
  partitioner.update_partition_pick_block(0, block_start, block_size); // modified by Jiang

  float gs[num_dims];
  for (int i = 0; i < num_dims; i ++) {
    bsz[i] = block_size[i];
    if (bsz[i] == _domain_size[i]) {
      gs[i] = 0;
      gsz[i] = bsz[i];
    } else {
      gs[i] = 1;
      gsz[i] = bsz[i] + 1;
    }
  }

  int left[num_dims];
  int min = 100000, min_dim = -1;
  for (int i = 0; i < num_dims; i ++) {  // find the shortest dimension
    left[i] = _domain_size[i] - gsz[i] ;
    if (left[i] < min && left[i] != 0) {
      min = left[i];
      min_dim = i;
    }
  }

  float ratio[num_dims];
  for (int i = 0; i < num_dims; i ++)  // ratio between each dimension and the shortest dimension
    ratio[i] = (float)left[i] / (float)left[min_dim];

  while (true) {  // compute max data with ghost layer (data duplication)
    bool flag = true;
    for (int i = 0; i < num_dims; i ++) {
      gsz[i] = std::min(bsz[i] + (int)std::ceil(gs[i]), _domain_size[i]);  // TODO 
      if ((float)bsz[i] + gs[i] < _domain_size[i]) flag = false;
    }
    
    if (flag) break;
    
    int64_t mem = 0, count = 1;
    for (int i = 0; i < num_dims; i ++) count *= gsz[i];
    count *= _dataconf.variables_size() * _dataconf.runs_size();
    mem += count * sizeof(float);
    
    if (mem * factor * (space_only() ? _domain_size[3] : 1) > mem_limit) break;

    for (int i = 0; i < num_dims; i ++) gs[i] += ratio[i];
  }
}

void CPTApp::probe_bs(int *bsz) const
{
  RegularPartitioner partitioner(
      _comm_world, _num_dims, _dataconf.runs_size(), _dataconf.variables_size());
  partitioner.set_domain_size(_domain_size);
  partitioner.set_num_blocks(_num_blocks);
  int block_start[4], block_size[4];
  partitioner.update_partition_pick_block(0, block_start, block_size); // modified by Jiang

  for (int i = 0; i < _num_dims; i ++)
    bsz[i] = block_size[i];
}

void CPTApp::probe_elastic_ghosts( 
    int *bsz, int *ghosts
  ) const // new refined by Jiang
{
  const int64_t mem_limit = (int64_t)_appconf.block_mem_limit() * 1024 * 1024 / (_num_blocks / comm_world_size());

  int num_dims = space_only() ? 3 : _num_dims;

  RegularPartitioner partitioner(
      _comm_world, num_dims, _dataconf.runs_size(), _dataconf.variables_size());
  partitioner.set_domain_size(_domain_size);
  partitioner.set_num_blocks(_num_blocks);
  int block_start[4], block_size[4];
  partitioner.update_partition_pick_block(0, block_start, block_size); // modified by Jiang

  int gsz[num_dims]; float gs[num_dims];
  for (int i = 0; i < num_dims; i ++) {
    bsz[i] = block_size[i];
    if (bsz[i] == _domain_size[i]) {
      gs[i] = 0;
      ghosts[i] = 0;
      gsz[i] = bsz[i];
    } else {
      gs[i] = 1;
      ghosts[i] = 1;
      gsz[i] = bsz[i] + 2;
    }
  }

  if (_appconf.block_mem_limit() < 0) {
    for (int i = 0; i < num_dims; i ++) ghosts[i] = 10000;
    return;
  } 

  int left[num_dims];
  int min = 100000, max = -1, min_dim = -1, max_dim = 5;
  for (int i = 0; i < num_dims; i ++) {  // find the shortest dimension
    left[i] = _domain_size[i] - gsz[i] ;
    if (left[i] < min && left[i] != 0) {
      min = left[i];
      min_dim = i;
    }

    if (left[i] > max && left[i] != 0) {
      max = left[i];
      max_dim = i;
    }
  }

  float ratio[num_dims];
  for (int i = 0; i < num_dims; i ++)  // ratio between each dimension and the shortest dimension
    ratio[i] = (float)left[i] / (float)left[max_dim];

  while (true) {  // compute max data with ghost layer (data duplication)
    bool flag = true;
    for (int i = 0; i < num_dims; i ++) {
      ghosts[i] = (int)std::ceil(gs[i]);
      gsz[i] = std::min(bsz[i] + ghosts[i] * 2, _domain_size[i]);  // TODO 
      if (gsz[i] < _domain_size[i]) flag = false;
    }
    
    if (flag) break;
    
    int64_t mem = 0, count = 1;
    for (int i = 0; i < num_dims; i ++) count *= gsz[i];
    count *= _dataconf.variables_size() * _dataconf.runs_size();
    mem += count * sizeof(float);
    
    if (mem * (space_only() ? _domain_size[3] : 1) > mem_limit) break;

    for (int i = 0; i < num_dims; i ++) gs[i] += ratio[i];
  }
}

void CPTApp::regular_probe_elastic_ghost_size(int *gs) const 
{
  const int64_t mem_limit = (int64_t)_appconf.block_mem_limit() * 1024 * 1024;

  RegularPartitioner partitioner(
      _comm_world, _dataconf.dimensions_size(), _dataconf.runs_size(), _dataconf.variables_size());
  partitioner.set_domain_size(_domain_size);
  partitioner.set_num_blocks(_num_blocks);
  partitioner.update_partition();

  for (int i = 0; i < num_dims(); i ++) gs[i] = 1;

  int j = 0, d = 0;
  while (true) {
    partitioner.set_ghost(gs);
    partitioner.update_ghosts(false);
    int64_t mem = partitioner.tot_block_mem();

    if (mem > mem_limit) break;
    if (partitioner.is_max_bounds()) break; // if not wrap boundary, should add this.

    int d = j % num_dims();
    gs[d] ++;
    j ++;
  }
}

void CPTApp::init_mpi()
{
  MPI_Comm_rank(_comm_world, &_comm_world_rank);
  MPI_Comm_size(_comm_world, &_comm_world_size);
}

void CPTApp::deinit()
{
  MPI_Barrier(_comm_world);
  _timestamp_finish = MPI_Wtime(); 
  _timestamps.push_back(_timestamp_finish); // wait
  _timecategories.push_back(2);
  stat(); 

#if WRITE
  write_output_file();
#endif
}

void CPTApp::stat()
{
  if (comm_world_rank() == 0) {
    fprintf(stderr, "--------------------------------------\n"); 
    fprintf(stderr, "--------------- SUMMARY --------------\n"); 
    fprintf(stderr, "time_init=\t\t%f\n", _timestamp_init - _timestamp_start); 
    fprintf(stderr, "time_run=\t\t%f\n", _timestamp_finish - _timestamp_init); 
    fprintf(stderr, "time_total=\t\t%f\n", _timestamp_finish - _timestamp_start); 
    fprintf(stderr, "--------------------------------------\n"); 
  }
}

void CPTApp::write_output_file()
{
  //fprintf(stderr, "rank %d: %d intervals\n", comm_world_rank(), (int)_timecategories.size());
  std::stringstream stream;  

  if (_appconf.block_mem_limit()>-1){
  stream << std::to_string(_appconf.block_mem_limit())<<"/gantt_"<<std::to_string(comm_world_size())<<"_"<< comm_world_rank() << ".csv"; 
  }else{
    stream << "unlim/gantt_"<<std::to_string(comm_world_size())<<"_"<< comm_world_rank() << ".csv";
  }
  std::ofstream ofile;  
  ofile.open(stream.str().c_str());
 
  int size = _timestamps.size();

  for (int i = 0; i < size-1; i ++) {
    float t = _timestamps[i+1] - _timestamps[i];
    if (i == size - 2) 
      ofile << _timecategories[i] << ":" << t;
    else 
      ofile << _timecategories[i] << ":" << t << ",";
  }

  ofile.close();

  fprintf(stderr, "rank %d write gantt finished\n", comm_world_rank());

  if (comm_world_rank() == 0) {
    std::ofstream ofile2;  
    stream.str("");
    if (_appconf.block_mem_limit()>-1){
    stream << std::to_string(_appconf.block_mem_limit())<<"/balance_"<<std::to_string(comm_world_size())<<".csv";
    }else{
      stream << "unlim/balance_"<<std::to_string(comm_world_size())<<".csv";
    }
    ofile2.open(stream.str().c_str());
 
    float prev = _timestamps[0]; int n = 0;
    //fprintf(stderr, "_timestamps.size=%d, _timecategories.size=%d, _balance.size=%d\n", 
    //  _timestamps.size(), _timecategories.size(), _balance.size());
    for (int i = 0; i < _timecategories.size(); i ++) {
      if (_timecategories[i] == 2) {
        ofile2 << _timestamps[i+1] - _timestamps[0] << ":" << _balance[n] << ":"<< _all_time_round_num[n]-1<< ",";
        n ++;
        prev = _timestamps[i+1];
      }
    }
    ofile2.close();

    fprintf(stderr, "rank %d write balance finished\n", comm_world_rank());


    stream.str("");
    if (_appconf.block_mem_limit()>-1){
    stream << std::to_string(_appconf.block_mem_limit())<<"/round_balance_"<<std::to_string(comm_world_size())<<".csv";
    }else{
      stream << "unlim/round_balance_"<<std::to_string(comm_world_size())<<".csv";
    }
    ofile2.open(stream.str().c_str());

    for (int i = 0; i < _round_balance.size(); i ++) {
        ofile2 << i << ":" << _round_balance[i] << ",";
    }

    ofile2.close();

    fprintf(stderr, "rank %d write round_balance finished\n", comm_world_rank());
  }
}

bool CPTApp::parse_arguments(int argc, char **argv)
{
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <config.xml>\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  if (!parse_config(argv[1])) { // default
    fprintf(stderr, "cannot parse %s\n", argv[1]);
    MPI_Abort(_comm_world, 0);
  }

  return true;
}

void CPTApp::bcast_config()
{
  int config_string_size=0, dataset_string_size=0; 
  std::string config_string, dataset_string; 
  if (comm_world_rank() == 0) {
    //parse_commandline(*argc, *argv); 
    _appconf.SerializeToString(&config_string); 
    config_string_size = config_string.size(); 
    _dataconf.SerializeToString(&dataset_string); 
    dataset_string_size = dataset_string.size(); 
  } 
  MPI_Bcast(&config_string_size, 1, MPI_INT, 0, _comm_world); 
  MPI_Bcast(&dataset_string_size, 1, MPI_INT, 0, _comm_world); 
  config_string.resize(config_string_size); 
  dataset_string.resize(dataset_string_size); 
  MPI_Bcast((char*)config_string.data(), config_string_size, MPI_CHAR, 0, _comm_world); 
  MPI_Bcast((char*)dataset_string.data(), dataset_string_size, MPI_CHAR, 0, _comm_world); 
  if (comm_world_rank() != 0) {
    _appconf.ParseFromString(config_string); 
    _dataconf.ParseFromString(dataset_string); 
  }

  _num_dims = _dataconf.dimensions_size();
  _nb_per_proc = _appconf.nb_per_proc();
  //if (is_baseline()) 
  { // TODO
    if (_nb_per_proc >= 1) _num_blocks = _nb_per_proc * comm_world_size();
    else _num_blocks = _appconf.num_blocks(); 
  } 
  //else
  //  _num_blocks = comm_world_size();
  _num_threads = _appconf.num_threads();

  for (int i = 0; i < num_dims(); i ++) 
    _domain_size[i] = _dataconf.dimensions(i).size();
}

bool CPTApp::parse_config(const std::string& filename)
{
  TiXmlDocument doc(filename.c_str()); 
  if (!doc.LoadFile()) return false; 

  TiXmlHandle hDoc(&doc); 
  TiXmlElement *pElem; 

  // app config
  { 
    pElem = hDoc.FirstChild("config").Element(); 
    assert(pElem != NULL); 
    
    // number of blocks in total (TODO move to number of blocks per proc.)
    int _n_blocks = _appconf.num_blocks();
    pElem->QueryIntAttribute("num_blocks", &_n_blocks);
    _appconf.set_num_blocks(_n_blocks);

    int _nb_proc = _appconf.nb_per_proc();
    pElem->QueryIntAttribute("nb_per_proc", &_nb_proc);
    _appconf.set_nb_per_proc(_nb_proc);

    int _epoch_size = _appconf.epoch_size(); 
    pElem->QueryIntAttribute("epoch_size", &_epoch_size); 
    _appconf.set_epoch_size(_epoch_size); 

    int _ghost_size = _appconf.ghost_size(); 
    pElem->QueryIntAttribute("ghost_size", &_ghost_size); 
    _appconf.set_ghost_size(_ghost_size); 

    int _num_threads = _appconf.num_threads();
    pElem->QueryIntAttribute("num_threads", &_num_threads); 
    _appconf.set_num_threads(_num_threads); 

    //if (pElem->Attribute("elastic_ghost")) 
    //  _appconf.set_elastic_ghost(true); 

    if (pElem->Attribute("elastic_ghost")) {
      int _elastic_ghost = _appconf.elastic_ghost();
      pElem->QueryIntAttribute("elastic_ghost", &_elastic_ghost);
      if (_elastic_ghost) _appconf.set_elastic_ghost(true);
      else _appconf.set_elastic_ghost(false);
    }

    if (pElem->Attribute("kd_tree")) {
      int _kd_tree = _appconf.kd_tree();
      pElem->QueryIntAttribute("kd_tree", &_kd_tree);
      if (_kd_tree) _appconf.set_kd_tree(true);
      else _appconf.set_kd_tree(false);
    } 

    if (pElem->Attribute("prediction")) {
      int _prediction = _appconf.prediction();
      pElem->QueryIntAttribute("prediction", &_prediction);
      _appconf.set_prediction(_prediction);
    } 

    /*if (pElem->Attribute("baseline")) {
      int _baseline = _appconf.baseline();
      pElem->QueryIntAttribute("baseline", &_baseline);
      if (_baseline) _appconf.set_baseline(true);
      else _appconf.set_baseline(false);
    }*/

    if (pElem->Attribute("space_only")) {
      int _space_only = _appconf.space_only();
      pElem->QueryIntAttribute("space_only", &_space_only);
      if (_space_only) _appconf.set_space_only(true);
      else _appconf.set_space_only(false);
    }
    
    int _block_mem_limit = _appconf.block_mem_limit();
    pElem->QueryIntAttribute("block_mem_limit", &_block_mem_limit);
    _appconf.set_block_mem_limit(_block_mem_limit);

    if (pElem->Attribute("debug")) 
      _appconf.set_debug(true); 
  }

  // data set config
  { 
    pElem = hDoc.FirstChild("dataset").Element(); 
    assert(pElem != NULL); 

    if (pElem->Attribute("name")) 
      _dataconf.set_name(pElem->Attribute("name")); 

    for (TiXmlElement *pElemDimension = pElem->FirstChild("dimension")->ToElement(); 
         pElemDimension != NULL; 
         pElemDimension = pElemDimension->NextSiblingElement("dimension"))
    {
      PBDataset::PBDimension *dimension = _dataconf.add_dimensions(); 

      if (pElemDimension->Attribute("name"))
        dimension->set_name(pElemDimension->Attribute("name")); 

      if (pElemDimension->Attribute("storage_order")) {
        std::string order = pElemDimension->Attribute("storage_order"); 
        if (order == "1") 
          dimension->set_storage_order(PBDataset::FIRST); 
        else if (order == "2") 
          dimension->set_storage_order(PBDataset::SECOND); 
        else if (order == "3") 
          dimension->set_storage_order(PBDataset::THIRD); 
        else if (order == "4") 
          dimension->set_storage_order(PBDataset::FOURTH); 
        else if (order == "FILE")
          dimension->set_storage_order(PBDataset::FILE); 
        else assert(false); 
      }

      if (pElemDimension->Attribute("size")) {
        int size; 
        pElemDimension->QueryIntAttribute("size", &size); 
        dimension->set_size(size); 
      }
    }

    for (TiXmlElement *pElemVariable = pElem->FirstChild("variable")->ToElement(); 
         pElemVariable != NULL; 
         pElemVariable = pElemVariable->NextSiblingElement("variable"))
    {
      PBDataset::PBVariable *variable = _dataconf.add_variables();

      if (pElemVariable->Attribute("name"))
        variable->set_name(pElemVariable->Attribute("name")); 

      if (pElemVariable->Attribute("data_type")) {
        std::string type = pElemVariable->Attribute("data_type"); 
        if (type == "FLOAT32") 
          variable->set_data_type(PBDataset::FLOAT32);
        else assert(false); 
      }

      if (pElemVariable->Attribute("big_endian")) {
        variable->set_big_endian(true); 
      }
    }

    for (TiXmlElement *pElemRun = pElem->FirstChild("run")->ToElement(); 
         pElemRun != NULL; 
         pElemRun = pElemRun->NextSiblingElement("run"))
    {
      PBDataset::PBRun *run = _dataconf.add_runs();

      if (pElemRun->Attribute("run_id"))
        run->set_run_id(pElemRun->Attribute("run_id")); 

      if (pElemRun->Attribute("directory")) 
        run->set_directory(pElemRun->Attribute("directory")); 

      if (pElemRun->Attribute("pattern")) 
        run->set_pattern(pElemRun->Attribute("pattern")); 
      
      std::string pattern = run->directory(); 
      pattern.append("/"); 
      pattern.append(run->pattern()); 

      glob_t results; 
      glob(pattern.c_str(), 0, NULL, &results); 
      for (int i=0; i<results.gl_pathc; i++) 
        run->add_filenames(results.gl_pathv[i]); 
      globfree(&results);

      if (results.gl_pathc == 0) {
        fprintf(stderr, "FATAL: cannot find data file %s\n", pattern.c_str());
        MPI_Abort(MPI_COMM_WORLD, 0);
      }
    }

    ////////
    PBDataset::PBRun *run0 = _dataconf.mutable_runs(0);
    int total_timesteps=0; 
    if (_dataconf.dimensions_size()==4){ //  && run0->file_names_size()>1) { // probe timesteps per file; // TODO
      for (int i=0; i<run0->filenames_size(); i++) {
        int ncid; 
        int varid; 
        int dimids[4]; 
        MPI_Offset dimlen; 
        
        PNC_SAFE_CALL( ncmpi_open(MPI_COMM_SELF, run0->filenames(i).c_str(), NC_NOWRITE, MPI_INFO_NULL, &ncid) ); 
        PNC_SAFE_CALL( ncmpi_inq_varid(ncid, _dataconf.variables(0).name().c_str(), &varid) ); 
        PNC_SAFE_CALL( ncmpi_inq_vardimid(ncid, varid, dimids) ); 
        PNC_SAFE_CALL( ncmpi_inq_dimlen(ncid, dimids[0], &dimlen) ); 

        // printf("%s, %d\n", _dataconf.file_names(i).c_str(), dimlen); 
        _dataconf.add_timesteps_per_file(dimlen);
        total_timesteps += dimlen; 

        PNC_SAFE_CALL( ncmpi_close(ncid) );
      }
      _dataconf.mutable_dimensions(3)->set_size(total_timesteps);
    }

    // fprintf(stderr, "num_data_files=%d, num_time_steps=%d\n", run0->file_names_size(), total_timesteps); 
  } 
  
  return true;
}

bool CPTApp::is_ptinblock(const Block& b, const float* pt)
{
  for (int i = 0; i < num_dims(); i ++) {
    if (b.ub[i] == domain_size()[i]-1) {
      if (pt[i] < b.lb[i] || pt[i] > b.ub[i]) return false;
    } else {
      if (pt[i] < b.lb[i] || pt[i] >= b.ub[i]+1) return false;
    }
  }

  return true;
}

int CPTApp::idx2id(const int idx[4]) const
{
  return idx[0] + _domain_size[0] * (idx[1] + _domain_size[1] * idx[2]); 
}

int CPTApp::max_num(const std::vector<int>& nums)
{
  int max = -1;
  for (int i = 0; i < nums.size(); i ++) {
    if (nums[i] > max) max = nums[i];
  }

  return max;
}

int CPTApp::min_num(const std::vector<int>& nums)
{
  int min = INT_MAX;
  for (int i = 0; i < nums.size(); i ++) {
    if (nums[i] < min) min = nums[i];
  }

  return min;
}

float CPTApp::avg_num(const std::vector<int>& nums)
{
 int64_t total = 0;
  for (int i = 0; i < nums.size(); i ++)
    total += nums[i];
  double avg = total/(double)nums.size();  

  return (float)avg;
}

float CPTApp::var_num(const std::vector<int>& nums)
{
  float avg = avg_num(nums), var = 0.;
  for (int i = 0; i < nums.size(); i ++)
    var += pow((float)nums[i]-avg, 2);

  return sqrt(var);
}

int CPTApp::sum_num(const std::vector<int>& nums)
{
  int total = 0;
  for (int i = 0; i < nums.size(); i ++)
    total += nums[i];

  return total;
}
