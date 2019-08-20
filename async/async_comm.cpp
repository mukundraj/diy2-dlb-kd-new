#include "async_comm.h"
#include <list>
#include <cassert>
#include <limits.h>

// from diy2-sfm

//static const int mpi_tag_shutdown = 3333; 

AsyncComm_MY::AsyncComm_MY(MPI_Comm comm, int max_nonblocking_requests) :
  _comm_world(comm),
  _max_nonblocking_requests(max_nonblocking_requests),
  _nonblocking_requests(max_nonblocking_requests, MPI_REQUEST_NULL),
  _nonblocking_buffers(max_nonblocking_requests),
  _num_blocked_sends(0),
  _num_bytes_sent(0)
{
  MPI_Comm_size(comm, &_np);
  MPI_Comm_rank(comm, &_rank);

  //for (int i=0; i<_max_nonblocking_requests; i++)
  //  _nonblocking_available.insert(i);
}

AsyncComm_MY::~AsyncComm_MY()
{
  // TODO
}

void AsyncComm_MY::reset_state()
{
  int all_finished = false, finished = false; 
  while (!all_finished) {
    recv_incoming_request();  
    finished = !has_nonblocking_requests(); 
    
    MPI_Allreduce(
        &finished, 
        &all_finished, 
        1, 
        MPI_INT, 
        MPI_LAND, 
        _comm_world
    ); 
    MPI_Barrier(_comm_world); 
    // printf("rank=%d\tfinished=%d\tall_finished=%d\n", comm_world_rank(), finished, all_finished); 
  }

  _nonblocking_requests.clear();
  _nonblocking_requests.resize(_max_nonblocking_requests);
  for (int i = 0; i <_max_nonblocking_requests; i ++)
    _nonblocking_requests[i] = MPI_REQUEST_NULL;
  _nonblocking_buffers.clear();
  _nonblocking_buffers.resize(_max_nonblocking_requests);

  MPI_Barrier(_comm_world);
}

void AsyncComm_MY::recv_incoming_request()
{
  MPI_Request request; 
  MPI_Status status; 
  int flag = 0; 

  MPI_Iprobe(
      MPI_ANY_SOURCE, 
      MPI_ANY_TAG, 
      _comm_world, 
      &flag, 
      &status
  ); 

  if (flag) {
    int length; 
    std::string buffer; 

    MPI_Get_count(&status, MPI_CHAR, &length); 
    buffer.resize(length); 

    MPI_Irecv(
      (void*)buffer.data(), 
      length, 
      MPI_CHAR, 
      status.MPI_SOURCE, 
      status.MPI_TAG, 
      _comm_world, 
      &request
    ); 
    MPI_Wait(&request, &status); 
  }
}

bool AsyncComm_MY::has_nonblocking_requests() 
{
  int flag; 
  for (std::vector<MPI_Request>::iterator it = _nonblocking_requests.begin(); 
       it != _nonblocking_requests.end(); 
       it ++) {
    MPI_Test(&(*it), &flag, MPI_STATUS_IGNORE); 
    if (*it != MPI_REQUEST_NULL) return true; 
  }
  return false; 
}

void AsyncComm_MY::wait_nonblocking_requests()
{
  for (std::vector<MPI_Request>::iterator it = _nonblocking_requests.begin(); 
       it != _nonblocking_requests.end(); 
       it ++) {
    // MPI_Request r = *it; 
    if (*it != MPI_REQUEST_NULL) { 
      int flag; 
      MPI_Test(&(*it), &flag, MPI_STATUS_IGNORE); 
      if (flag) {
        // fprintf(stderr, "[%d] waiting for non-blocking requests...\n", comm_world_rank()); 
        MPI_Wait(&(*it), MPI_STATUS_IGNORE); 
        *it = MPI_REQUEST_NULL; 
      }
    }
  }
}

bool AsyncComm_MY::ready_to_send()
{
  return true; // TODO
}

bool AsyncComm_MY::isend(int dst_rank, int tag, diy::MemoryBuffer& bb)  
{
  if (bb.empty()) return false;

  int request = get_nonblocking_request(); 
  if (request == -1) {
    _num_blocked_sends ++; 
    return false; 
  }

  diy::MemoryBuffer &buffer = _nonblocking_buffers[request]; 
  buffer.swap(bb); 

  MPI_Isend(
    (void*)buffer.buffer.data(), 
    buffer.size(), 
    MPI_CHAR,
    dst_rank, 
    tag, 
    _comm_world, 
    &_nonblocking_requests[request]
  );

  // fprintf(stdout, "[rank=%d] sent %lu bytes\n", rank(), buffer.size());

  _num_bytes_sent += buffer.size();
  return true; 
}

bool AsyncComm_MY::iprobe(int &src_rank, int tag, diy::MemoryBuffer& bb)
{
  MPI_Request request; 
  MPI_Status status; 
  int flag = 0; 

  MPI_Iprobe(
    MPI_ANY_SOURCE, 
    tag, 
    _comm_world, 
    &flag, 
    &status
  ); 

  if (flag) { // has message 
    int length; 

    src_rank = status.MPI_SOURCE;

    MPI_Get_count(&status, MPI_CHAR, &length); 
    bb.clear();
    bb.buffer.resize(length); 
    assert(length != 0);

    MPI_Irecv(
      (void*)bb.buffer.data(), 
      length, 
      MPI_CHAR, 
      status.MPI_SOURCE, 
      status.MPI_TAG, 
      _comm_world, 
      &request
    ); 
    MPI_Wait(&request, &status); 
  
    // fprintf(stdout, "[rank=%d] recv %lu bytes\n", rank(), bb.size());

    return true; 
  } else 
    return false; 
}

int AsyncComm_MY::get_nonblocking_request()
{
  clear_finished_nonblocking_request(); 

  for (int i = 0; i < _max_nonblocking_requests; i ++) 
    if (_nonblocking_requests[i] == MPI_REQUEST_NULL)
      return i; 
  return -1;
}

void AsyncComm_MY::clear_finished_nonblocking_request()
{
  for (int i = 0; i < _max_nonblocking_requests; i ++) {
    if (_nonblocking_requests[i] != MPI_REQUEST_NULL) {
      int flag; // finished 
      MPI_Test(
        &_nonblocking_requests[i], 
        &flag, 
        MPI_STATUS_IGNORE
      );

      if (flag) {
        MPI_Wait(
          &_nonblocking_requests[i], 
          MPI_STATUS_IGNORE
        );
        _nonblocking_buffers[i].clear(); 
      }
    }
  }
}
