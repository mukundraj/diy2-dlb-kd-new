#ifndef _CPTAPP_ASYNC_H
#define _CPTAPP_ASYNC_H

#include "app.h"
#include "async/async_master.hpp"

#define NUM_ROUNDS 40 // 30 for nek; 20 for isabel; 30 for geos5

class CPTApp_Async : public CPTApp {
  friend class VerifyBlockAsync;
  friend class TraceBlockAsync;
  friend class TraceBlockAsync_Interval;

public:
  CPTApp_Async();
  ~CPTApp_Async();

  void init(int argc, char **argv);
  void exec();

  AsyncMaster<Message>& asyncMaster() { return *_async_master; }

public: // interface 
  virtual void initialize_particles(Block&,
      std::vector<Particle>&
    ) = 0;
  virtual void trace_particles(Block&, 
      std::vector<Particle>&, // particles to trace
      std::map<int, std::vector<Particle> >&, // unfinished particles
      std::map<int, std::vector<Particle> >& // finished particles
    ) = 0;
/*
  virtual void trace_particles_kdtree(Block&, 
      std::vector<Particle>&, // particles to trace
      std::map<int, std::vector<Particle> >&, // unfinished particles
      std::map<int, std::vector<Particle> >& // finished particles
    ) = 0;
*/
  virtual void trace_particles_kdtree(Block&, 
      std::vector<Particle>&, // particles to trace
      std::map<int, std::vector<Particle> >&, // unfinished particles
      std::map<int, int>& // number of finished particles
    ) = 0;

protected:
  void stat();
  mutable double _time_work, _time_kdtree, _time_exchange, _time_prediction;
  mutable double _time_trace, _time_iexchange;

protected:
  void add_workload();
  int _total_steps, _round_steps; 

  uint64_t _max_workload;
  int _num_particles;

  int _num_local_initialized, _num_local_finished;

private:
  mutable AsyncMaster<Message> *_async_master;

};

#endif
