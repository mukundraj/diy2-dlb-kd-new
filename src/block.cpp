#include "block.h"
#include "bil/bil.h"
#include <vector>
#include <string>
#include <mpi.h>
#include <algorithm>
#include <cassert>
#include <boost/thread/mutex.hpp>

static boost::mutex bil_mutex;

static void mybil_add_block( // BIL_Ext_Add_block_nc(
    int num_dims, 
    const int *block_starts, 
    const int *block_sizes, 
    const std::vector<std::string> &filenames, 
    const std::vector<int> timesteps_per_file,
    const std::string& varname,
    float *buf)
{
  assert(num_dims==4); 

  int t0=block_starts[0], 
      ts=block_sizes[0]; 

  int T0=0; 
  int i; 
  for (i=0; i<timesteps_per_file.size(); i++) {
    if (T0+timesteps_per_file[i] > t0) break; 
    else T0+=timesteps_per_file[i]; 
  }

  while (1) {
    int lt0 = t0-T0; 
    int lts = std::min(timesteps_per_file[i], ts); 

#if 1 // WITH_BIL
    int starts[4] = {lt0, block_starts[1], block_starts[2], block_starts[3]}, 
        sizes[4]  = {lts, block_sizes[1],  block_sizes[2],  block_sizes[3]};
    BIL_Add_block_nc(num_dims, starts, sizes, filenames[i].c_str(), varname.c_str(), (void**)&buf);
#else // raw netcdf
    int ncid; 
    int varid; 
    size_t starts[4] = {lt0, block_starts[1], block_starts[2], block_starts[3]}, 
           sizes[4]  = {lts, block_sizes[1],  block_sizes[2],  block_sizes[3]}; 
    
    NC_SAFE_CALL( nc_open(filenames[i].c_str(), NC_NOWRITE, &ncid) ); 
    NC_SAFE_CALL( nc_inq_varid(ncid, varname, &varid) ); 
    NC_SAFE_CALL( nc_get_vara_float(ncid, varid, starts, sizes, (float*)buf) ); 
    NC_SAFE_CALL( nc_close(ncid) ); 
#endif
    buf += sizes[0]*sizes[1]*sizes[2]*sizes[3]; 

    ts -= lts; 
    if (ts==0) break; 
    else i++; 
  }
}

Block::Block() :
  num_particles_initialized(-1), // very beginning
  num_particles_finished(0),
  num_particles_traced(0),
  num_particles_kd(0),
  //num_particles(0),
  initialized(false), 
  finished(false),
  workload(0),
  after_kdtree(false)
{
  memset(lb, 0, sizeof(int)*4);
  memset(ub, 0, sizeof(int)*4);
  memset(glb, 0, sizeof(int)*4);
  memset(gub, 0, sizeof(int)*4);
  memset(lload, 0, sizeof(int)*4);
  memset(uload, 0, sizeof(int)*4);

  //for (int i = 0; i < 8; i ++) 
  //  neighbors[i] = -1;
}
  
Block::~Block() {
  for (int i = 0; i < vars.size(); i ++) 
    free(vars[i]);
}

void Block::bil_add_block_2D(const PBDataset& d)
{
  const int starts[2] = {lload[1], lload[0]}, 
            sizes[2]  = {uload[1]-lload[1]+1, uload[0]-lload[0]+1};
  const int arraySize = sizes[0] * sizes[1];

  for (int r=0; r<d.runs_size(); r++) {
    const PBDataset::PBRun& run = d.runs(r);
    const std::string& filename = run.filenames(0); 
    for (int v=0; v<d.variables_size(); v++) {
      boost::mutex::scoped_lock lock(bil_mutex);
      const std::string& varname = d.variables(v).name(); 
      float *p = (float*)malloc(sizeof(float)*arraySize);
      vars.push_back(p);
      BIL_Add_block_nc(2, starts, sizes, filename.c_str(), varname.c_str(), (void**)&p);
    }
  }
}
  
void Block::bil_add_block_3D(const PBDataset& d)
{
  //const int starts[3] = {glb[2], glb[1], glb[0]}, 
  //          sizes[3]  = {gub[2]-glb[2]+1, gub[1]-glb[1]+1, gub[0]-glb[0]+1};
  const int starts[3] = {lload[2], lload[1], lload[0]}, 
            sizes[3]  = {uload[2]-lload[2]+1, uload[1]-lload[1]+1, uload[0]-lload[0]+1};
  const int arraySize = sizes[0] * sizes[1] * sizes[2];

  for (int r=0; r<d.runs_size(); r++) {
    const PBDataset::PBRun& run = d.runs(r);
    const std::string& filename = run.filenames(0); 
    for (int v=0; v<d.variables_size(); v++) {
      boost::mutex::scoped_lock lock(bil_mutex);
      const std::string& varname = d.variables(v).name(); 
      float *p = (float*)malloc(sizeof(float)*arraySize);
      vars.push_back(p);
      BIL_Add_block_nc(3, starts, sizes, filename.c_str(), varname.c_str(), (void**)&p);
    }
  }
}

void Block::bil_add_block_4D(const PBDataset& d)
{
  //const int starts[4] = {glb[3], glb[2], glb[1], glb[0]}, 
  //          sizes[4]  = {gub[3]-glb[3]+1, gub[2]-glb[2]+1, gub[1]-glb[1]+1, gub[0]-glb[0]+1};
  const int starts[4] = {lload[3], lload[2], lload[1], lload[0]}, 
            sizes[4]  = {uload[3]-lload[3]+1, uload[2]-lload[2]+1, uload[1]-lload[1]+1, uload[0]-lload[0]+1};
  const int arraySize = sizes[0] * sizes[1] * sizes[2] * sizes[3]; 

  for (int r=0; r<d.runs_size(); r++) {
    const PBDataset::PBRun& run = d.runs(r);
    std::vector<std::string> filenames; 
    std::vector<int> timesteps_per_file; 

    for (int i=0; i<run.filenames_size(); i++) {
      filenames.push_back(run.filenames(i)); 
      timesteps_per_file.push_back(d.timesteps_per_file(i)); 
    }

    for (int v=0; v<d.variables_size(); v++) {
      boost::mutex::scoped_lock lock(bil_mutex);
      const std::string& varname = d.variables(v).name(); 
      float *p = (float*)malloc(sizeof(float)*arraySize);
      vars.push_back(p);
      mybil_add_block(4, starts, sizes, filenames, timesteps_per_file, varname, p);
    }
  }
}

void Block::add_block_nek(const PBDataset& d)
{
  std::vector<std::string> varnames;
  varnames.push_back("U"); varnames.push_back("V"); varnames.push_back("W");
#if 1
  //MPI_Barrier(MPI_COMM_WORLD);
  //double start = MPI_Wtime();
  std::vector<std::string> filenames;
  filenames.push_back("/projects/SDAV/jiangz/nek5000/U.nc"); 
  filenames.push_back("/projects/SDAV/jiangz/nek5000/V.nc"); 
  filenames.push_back("/projects/SDAV/jiangz/nek5000/W.nc");

  const int starts[3] = {lload[2], lload[1], lload[0]}, 
            sizes[3]  = {uload[2]-lload[2]+1, uload[1]-lload[1]+1, uload[0]-lload[0]+1};
  const int arraySize = sizes[0] * sizes[1] * sizes[2];

  for (int i = 0; i < 3; i ++) {
    boost::mutex::scoped_lock lock(bil_mutex);
    const std::string& filename = filenames[i]; 
    const std::string& varname = varnames[i]; 
    float *p = (float*)malloc(sizeof(float)*arraySize);
    vars.push_back(p);
    BIL_Add_block_nc(3, starts, sizes, filename.c_str(), varname.c_str(), (void**)&p);
  }

  //MPI_Barrier(MPI_COMM_WORLD);
  //double end = MPI_Wtime();
  //iotime = end - start;
#else
  // TODO
#endif
}
  
void Block::show_geometry(const diy::Master::ProxyWithLink &cp, void*)
{
  fprintf(stderr, "rank=%d, gid=%d, lb={%d, %d, %d, %d}, ub={%d, %d, %d, %d}, glb={%d, %d, %d, %d}, gub={%d, %d, %d, %d}, lload={%d, %d, %d, %d}, uload={%d, %d, %d, %d}\n", 
      cp.master()->communicator().rank(),
      gid,
      lb[0], lb[1], lb[2], lb[3], 
      ub[0], ub[1], ub[2], ub[3],
      glb[0], glb[1], glb[2], glb[3], 
      gub[0], gub[1], gub[2], gub[3],
      lload[0], lload[1], lload[2], lload[3],
      uload[0], uload[1], uload[2], uload[3]);
}

void Block::print_info()
{
  fprintf(stderr, "lb={%d, %d, %d, %d}, ub={%d, %d, %d, %d}, glb={%d, %d, %d, %d}, gub={%d, %d, %d, %d}, lload={%d, %d, %d, %d}, uload={%d, %d, %d, %d}\n",
      lb[0], lb[1], lb[2], lb[3],
      ub[0], ub[1], ub[2], ub[3],
      glb[0], glb[1], glb[2], glb[3],
      gub[0], gub[1], gub[2], gub[3],
      lload[0], lload[1], lload[2], lload[3],
      uload[0], uload[1], uload[2], uload[3]);
}
  
void Block::get_ghost_st_sz(int ndims, int *gst, int *gsz)
{
  for (int i=0; i<ndims; i++) {
    gst[i] = glb[i]; 
    gsz[i] = gub[i] - glb[i] + 1;
  }
}

void Block::get_load_st_sz(int ndims, int *lst, int *lsz)
{
  for (int i=0; i<ndims; i++) {
    lst[i] = lload[i]; 
    lsz[i] = uload[i] - lload[i] + 1;
  }
}

void Block::get_ghost_load_st_sz(int ndims, int *gst, int *gsz, int *lst, int *lsz)
{
  for (int i=0; i<ndims; i++) {
    gst[i] = glb[i]; 
    gsz[i] = gub[i] - glb[i] + 1;
    lst[i] = lload[i]; 
    lsz[i] = uload[i] - lload[i] + 1;
  }
}

/*
int Block::neighbor_dir(int gid) const 
{
  for (int i=0; i<8; i++) 
    if (neighbors[i] == gid) return i;
  return -1;
}
*/
