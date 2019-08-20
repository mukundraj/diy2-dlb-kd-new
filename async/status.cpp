#include "status.h"

#define ASYNC 1

// from diy2-sfm 
// Francez’s algorithm for the detection of global termination

StatusExchange::StatusExchange(MPI_Comm comm, int num_local_blocks) :
  _comm(comm), 
  _num_local_blocks(num_local_blocks),
  _global_all_done(false),
  _local_all_done(false),
  _num_local_blocks_done(0),
  _num_ranks_done(0)
{
  MPI_Comm_size(comm, &_np);
  MPI_Comm_rank(comm, &_rank);
}

StatusExchange::~StatusExchange()
{
}

void StatusExchange::reset_state()
{
  _global_all_done = false;
  _local_all_done = false;
  _num_local_blocks_done = 0;
  _num_ranks_done = 0;
}

void StatusExchange::iexchange()
{
  if (_global_all_done) return;

  static int dummy = 0; // we do not care what number is sent
  const int tag_local_all_done = 3333, tag_global_all_done = 3334; // TODO
  MPI_Status status;
  MPI_Request request;
  int flag;
  
  MPI_Iprobe(MPI_ANY_SOURCE, tag_global_all_done, _comm, &flag, &status);
  if (flag) {
    int src_rank = status.MPI_SOURCE;
    MPI_Recv(&dummy, 1, MPI_INT, src_rank, tag_global_all_done, _comm, MPI_STATUS_IGNORE);

    if (_rank != 0 && _rank*2 < _np) 
#if ASYNC
      MPI_Isend(&dummy, 1, MPI_INT, _rank*2, tag_global_all_done, _comm, &request);
#else
      MPI_Send(&dummy, 1, MPI_INT, _rank*2, tag_global_all_done, _comm);
#endif
    if (_rank*2+1 < _np) 
#if ASYNC
      MPI_Isend(&dummy, 1, MPI_INT, _rank*2+1, tag_global_all_done, _comm, &request);
#else
      MPI_Send(&dummy, 1, MPI_INT, _rank*2+1, tag_global_all_done, _comm);
#endif
    
    _global_all_done = true; // game over~
  }

  // probing locall_all_done signals
  MPI_Iprobe(MPI_ANY_SOURCE, tag_local_all_done, _comm, &flag, &status);
  if (flag) {
    int src_rank = status.MPI_SOURCE;
    MPI_Recv(&dummy, 1, MPI_INT, src_rank, tag_local_all_done, _comm, MPI_STATUS_IGNORE);

    if (_rank != 0) {
#if ASYNC
      MPI_Isend(&dummy, 1, MPI_INT, _rank/2, tag_local_all_done, _comm, &request);
#else
      MPI_Send(&dummy, 1, MPI_INT, _rank/2, tag_local_all_done, _comm);
#endif
    } else { // i'm rank 0
      _num_ranks_done ++;
      if (_num_ranks_done == _np) { // global_all_done
        if (_np > 1)
#if ASYNC
          MPI_Isend(&dummy, 1, MPI_INT, 1, tag_global_all_done, _comm, &request);
#else
          MPI_Send(&dummy, 1, MPI_INT, 1, tag_global_all_done, _comm);
#endif
        _global_all_done = true; // game over~
      }
    }
  }

  if (!_local_all_done) {
    if (_num_local_blocks_done >= _num_local_blocks) { // _num_local_blocks_done == _num_local_blocks
      _local_all_done = true;
      // fprintf(stderr, "rank=%d, local all done!\n", _rank);
      if (_rank != 0) {
        MPI_Send(&_rank, 1, MPI_INT, _rank/2, tag_local_all_done, _comm);
      } else { // i'm rank 0
        _num_ranks_done ++;
        if (_num_ranks_done == _np) { // global_all_done
          if (_np > 1)
            MPI_Send(&dummy, 1, MPI_INT, 1, tag_global_all_done, _comm);
          _global_all_done = true; // game over~
        }
      }
    }
  }
}
