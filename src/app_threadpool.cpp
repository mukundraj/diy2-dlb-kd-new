#include <cassert>
#include <boost/foreach.hpp>
#include "app_threadpool.h"
#include "cons_kdtree.h"

CPTApp_ThreadPool::CPTApp_ThreadPool() : 
  CPTApp(),
  _time_work(0), 
  _time_exchange(0),
  _time_trace(0)
{
}

CPTApp_ThreadPool::~CPTApp_ThreadPool()
{
  if (_threadpool_master != NULL) delete _threadpool_master;
}

void CPTApp_ThreadPool::init(int argc, char **argv)
{
  CPTApp::init(argc, argv);
  _threadpool_master = new ThreadPoolMaster<Message, ProcStat>(comm_world(), *_assigner, gids().size());
}

void CPTApp_ThreadPool::exec()
{
  if (comm_world_rank() == 0) 
    fprintf(stderr, "running...\n");

  // if (!is_kd_tree()) {} // TODO

  // set cbs
  _threadpool_master->set_msg_cb(boost::bind(&CPTApp_ThreadPool::process_message, this, _1, _2));
  //_threadpool_master->set_pstat_cbs(
  //    boost::bind(&CPTApp_ThreadPool::get_pstat, this, _1), 
  //    boost::bind(&CPTApp_ThreadPool::proc_pstat, this, _1, _2));
  //_tpm->set_cycle_cb(boost::bind(&ParticleTracingApp_ThreadPool::balance_workload, this));

  BOOST_FOREACH(int gid, gids()) {
    Message m;
    enqueue(gid, gid, Message::MESSAGE_INIT, m);
  }

  // run threads
  _threadpool_master->run(num_threads()); 
}

void CPTApp_ThreadPool::stat()
{ 
  // print workload
  std::vector<uint64_t> blk_workloads(_num_blocks, 0);
  uint64_t proc_workloads = 0;
  uint64_t proc_num_particles = 0;
  BOOST_FOREACH(int gid, gids()) {
    Block &b = block(gid);
    blk_workloads[gid] = b.workload;
    proc_workloads += b.workload;
    proc_num_particles += b.num_particles_finished;
  }
 
  uint64_t all_num_particles = 0;
  MPI_Reduce(
      &proc_num_particles, 
      &all_num_particles, 
      1, MPI_UINT64_T, MPI_SUM, 0, _comm_world);

  uint64_t proc_num_bytes_sent = _threadpool_master->num_bytes_sent(), 
           all_num_bytes_sent = 0;
  MPI_Reduce(
      &proc_num_bytes_sent, 
      &all_num_bytes_sent, 
      1, MPI_UINT64_T, MPI_SUM, 0, _comm_world);

  std::vector<uint64_t> all_blk_workloads(_num_blocks, 0);
  MPI_Reduce(
      (uint64_t*)blk_workloads.data(), // sendbuf
      (uint64_t*)all_blk_workloads.data(), // recvbuf
      _num_blocks, // count
      MPI_UINT64_T, MPI_SUM, 
      0, _comm_world);
  
  std::vector<uint64_t> all_proc_workloads(comm_world_size());
  MPI_Gather(&proc_workloads, 1, MPI_UINT64_T,
      (uint64_t*)all_proc_workloads.data(), 1, MPI_UINT64_T, 0, _comm_world);
   
  if (comm_world_rank() == 0) {
    // for (int i=0; i<all_blk_workloads.size(); i++)
    //   fprintf(stderr, "gid=%d, workload=%lld\n", i, all_blk_workloads[i]);
    fprintf(stderr, "num_bytes_sent=%llu\n", all_num_bytes_sent);
    fprintf(stderr, "num_particles=%llu\n", all_num_particles);
    fprintf(stderr, "imba_blks=%.4f%%\n", imbalance(all_blk_workloads)*100);
    fprintf(stderr, "imba_proc=%.4f%%\n", imbalance(all_proc_workloads)*100);
  }
}

void CPTApp_ThreadPool::process_message(int gid, Message& m)
{
  Block &b = block(gid);
  assert(b.gid == gid);
  
  std::map<int, std::vector<Particle> > unfinished_particles, finished_particles;

  if (m.type == Message::MESSAGE_INIT) {
    // fprintf(stderr, "[%d] INIT\n", gid);
    Message outmsg;
    initialize_particles(b, outmsg.unfinished_particles); // TODO 

    //const int nkeys = b.idmap.size();
    //b.initialized_particle_counts.resize(nkeys, num_sto_runs());
    //b.finished_particle_counts.resize(nkeys, 0);
    //b.finished_particles.resize(nkeys);
    //for (int i=0; i<nkeys; i++) 
    //  b.finished_particles[i].resize(num_sto_runs());
    
    // for (int i=0; i<num_sto_runs()-1; i++)
    //  outmsg.unfinished_particles.insert(
    //      outmsg.unfinished_particles.end(), 
    //      outmsg.unfinished_particles.begin(), 
    //      outmsg.unfinished_particles.begin() + b.idmap.size());
    
    b.num_particles_initialized += outmsg.unfinished_particles.size();
    enqueue(gid, gid, Message::MESSAGE_TRACE, outmsg);
  } else if (m.type == Message::MESSAGE_TRACE) {
    // fprintf(stderr, "[%d] TRACE, %d\n", gid, m.unfinished_particles.size());
    const int count = m.unfinished_particles.size();
    trace_particles(b, m.unfinished_particles, unfinished_particles, finished_particles);
    
    // fprintf(stderr, "task_gpu=%d, n=%d, t=%.2f\n", gpu, (int)m.unfinished_particles.size(), t1-t0);
    
    for (std::map<int, std::vector<Particle> >::iterator it = unfinished_particles.begin(); it != unfinished_particles.end(); it ++) {
      Message outmsg;
      outmsg.unfinished_particles.swap(it->second);
      enqueue(gid, it->first, Message::MESSAGE_TRACE, outmsg);
      assert(gid != it->first);
    }
    for (std::map<int, std::vector<Particle> >::iterator it = finished_particles.begin(); it != finished_particles.end(); it ++) {
      Message outmsg;
      outmsg.finished_particles.swap(it->second);
      enqueue(gid, it->first, Message::MESSAGE_FINISH, outmsg);
    }
  } else if (m.type == Message::MESSAGE_FINISH) {
    b.num_particles_finished += m.finished_particles.size(); 
  }
    
  // check if queries are finished.  only one thread can change status
  // fprintf(stderr, "[%d] checking, init=%d, done=%d\n", gid, (int)b.num_particles_initialized, (int)b.num_particles_finished);
  if (!b.finished && b.num_particles_initialized == b.num_particles_finished && b.num_particles_finished != 0) {
    if (b.mutex_finish.try_lock()) {
      b.finished = true;
      // _finished_blocks.insert(b.gid);
      _threadpool_master->inc_local_blocks_done(); 
      fprintf(stderr, "[rank=%d] blk %d done, #particles=%d\n", 
          comm_world_rank(), b.gid, (int)b.num_particles_finished);
      b.mutex_finish.unlock();
    } else {
      b.mutex_finish.lock();
      b.mutex_finish.unlock();
    }
  }
}

void CPTApp_ThreadPool::enqueue(int src_gid, int dst_gid, int type, Message& m)
{
  const int limit = 64; //TODO //task_max_cpu();
  if (m.unfinished_particles.size() > limit) {
    const int count = m.unfinished_particles.size(); 
    const int nmsg = (count%limit == 0) ? (count/limit) : (count/limit+1);

    for (int i=0; i<nmsg; i++) {
      Message m1;
      int size; 
      if (i < nmsg-1) size = limit;
      else size = count - limit*i;

      m1.unfinished_particles.resize(size);
      memcpy(
          &m1.unfinished_particles[0],
          &m.unfinished_particles[i*limit],
          size * sizeof(Particle));
      enqueue(src_gid, dst_gid, type, m1);
    }
  } else {
    m.type = type;
    //m.src_workload = block(src_gid).workload;
    _threadpool_master->enqueue(src_gid, dst_gid, m);
    // fprintf(stderr, "real enqueue, src_gid=%d, dst_gid=%d, type=%d, size=%d, time=%f\n", 
    //     src_gid, dst_gid, type, m.unfinished_particles.size(), t1-t0);
  }
}
