#ifndef _STATUS_EXCHANGE_H
#define _STATUS_EXCHANGE_H

#include <mpi.h>

class StatusExchange {
public:
  StatusExchange(MPI_Comm comm, int num_local_blocks);
  ~StatusExchange();

  void iexchange(); // called from main thread
  
  bool global_all_done() const { return _global_all_done; }
  
  void inc_local_blocks_done() { _num_local_blocks_done ++; }
  void set_local_blocks_done(int n) { _num_local_blocks_done = n; }

  void reset_state();

private:
  bool _global_all_done, _local_all_done;
  int _num_local_blocks_done, _num_ranks_done;
  
  const int _num_local_blocks;
  MPI_Comm _comm;
  int _np, _rank;
  
};

#endif
