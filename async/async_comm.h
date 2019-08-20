#ifndef _ASYNC_COMM_MY_H
#define _ASYNC_COMM_MY_H

#include <mpi.h>
#include <vector>
#include <string>
#include <diy/serialization.hpp>
#include <stdint.h>

// from diy2-sfm

class AsyncComm_MY { // not thread-safe. 
public:
  AsyncComm_MY(MPI_Comm comm=MPI_COMM_WORLD, int max_nonblocking_requests=1024);
  ~AsyncComm_MY();

  int rank() const { return _rank; }
  int np() const { return _np; }

  bool ready_to_send(); // ready to send
  bool isend(int dst_rank, int tag, diy::MemoryBuffer&);
  bool iprobe(int &src_rak, int tag, diy::MemoryBuffer&);

  void reset_state();

  uint64_t num_bytes_sent() const { return _num_bytes_sent; }

private:
  int get_nonblocking_request(); 
  //void clear_finished_nonblocking_request(int num=-1);
  void clear_finished_nonblocking_request();

  void recv_incoming_request();
  bool has_nonblocking_requests();
  void wait_nonblocking_requests();

private:
  const int _max_nonblocking_requests;
  std::vector<MPI_Request> _nonblocking_requests; 
  std::vector<diy::MemoryBuffer> _nonblocking_buffers; 
  //std::set<int> _nonblocking_available, _nonblocking_occupied;

protected:
  int _np, _rank;
  MPI_Comm _comm_world;
  uint64_t _num_blocked_sends, _num_bytes_sent;

};

#endif
