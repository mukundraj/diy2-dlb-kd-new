#ifndef _CPTAPP_SYNC_H
#define _CPTAPP_SYNC_H

#include "app.h"

class CPTApp_Sync : public CPTApp {
  friend class TraceBlockSync;
  friend class VerifyBlockSync;

public:
  CPTApp_Sync();

  void init(int argc, char **argv);
  void exec();

public: // interface 
  virtual void initialize_particles(Block&,
      std::vector<Particle>&
    ) = 0;

  virtual void trace_particles_core(Block&, 
      std::vector<Particle>&, // particles to trace
      std::map<int, std::vector<Particle> >&, // unfinished particles
      std::map<int, std::vector<Particle> >& // finished particles
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

  virtual void trace_particles_kdtree_predict(Block&, int) = 0;

  virtual void trace_particles_video(Block&, 
      std::vector<Particle>&, 
      std::vector<Particle>&,
      std::vector<Particle>&
    ) = 0;

protected:
  void stat();
  mutable double _time_work, _time_kdtree, _time_exchange, _time_prediction;
  mutable double _time_trace;

  int _local_init, _local_done, _local_init_epoch, _local_done_epoch;

protected:
  void gather_store_cores(const Block& b);
  void gather_store_particles(const int count, const std::vector<Particle>& particles);
  void add_workload();
  int _total_steps, _round_steps, _pred_mismatch, _total_pred_mismatch; 

  uint64_t _max_workload;
  int _num_particles;

};

#endif
