#ifndef _APP_THREADPOOL_H
#define _APP_THREADPOOL_H

#include "app.h"
#include "async/threadpool_master.hpp"

class CPTApp_ThreadPool : public CPTApp {
  typedef float ProcStat; // load balancing 

public:
  CPTApp_ThreadPool();
  ~CPTApp_ThreadPool();
  
  static int mpi_thread_level() { return MPI_THREAD_FUNNELED; }

  void init(int argc, char **argv);
  void exec();
  // void deinit();

  ThreadPoolMaster<Message, ProcStat>& threadPoolMaster() { return *_threadpool_master; }

public: // interface 
  virtual void initialize_particles(Block&,
      std::vector<Particle>&
    ) = 0;

  virtual void trace_particles(Block&, 
      std::vector<Particle>&, // particles to trace
      std::map<int, std::vector<Particle> >&, // unfinished particles
      std::map<int, std::vector<Particle> >& // finished particles
    ) = 0;

  virtual void trace_particles_kdtree(Block&, 
      std::vector<Particle>&, // particles to trace
      std::map<int, std::vector<Particle> >&, // unfinished particles
      std::map<int, std::vector<Particle> >& // finished particles
    ) = 0;

  void inc_workload(Block& b, int n) { b.workload += n; }

protected: 
  // called from worker threads
  void process_message(int gid, Message& m);
  // called from the comm thread
  //bool process_message_instant(int gid, Message& m);

protected:
  void enqueue(int src_gid, int dst_gid, int type, Message& m);
  // void split_message(Message& m, std::vector<Message> &ms, int limit);

protected:
  void stat();
  mutable double _time_work, _time_trace, _time_exchange;

private:
  ThreadPoolMaster<Message, ProcStat> *_threadpool_master;
  // std::set<int> _finished_blocks; // key=gid

};

#endif
