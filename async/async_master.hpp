#ifndef _ASYNC_MASTER_HPP
#define _ASYNC_MASTER_HPP

#include "async_comm.h"
#include "comm.h"
#include "status.h"
#include <mpi.h>
#include <diy/assigner.hpp>
#include <vector>

//////////////////////// from diy2-sfm, but no thread (un)lock

#if 0
typedef AsyncComm_MY AsyncComm;
#else
typedef AsyncComm_HG AsyncComm;
#endif

template<class T>
class AsyncMaster { 
public:
  AsyncMaster(MPI_Comm comm, const diy::Assigner& assigner, int num_local_blocks) :
    // _comm(comm),
    _assigner(assigner),
    _status_exchange(comm, num_local_blocks),
    _async(comm, 1024)
  {}
  ~AsyncMaster() {}
  
  void inc_local_blocks_done() { _status_exchange.inc_local_blocks_done(); }
  void set_local_blocks_done(int n) { _status_exchange.set_local_blocks_done(n); }
  //void reset() { _inflight.clear(); _outflight.clear(); _async.reset_comm(); _status_exchange.reset_status(); }
  bool running() const { return !_status_exchange.global_all_done(); }

  void reset_state() { _inflight.clear(); _outflight.clear(); _async.reset_state(); _status_exchange.reset_state(); }
 
  void enqueue(
      int src_gid,
      int dst_gid, 
      T& x, 
      void (*save)(diy::BinaryBuffer&, const T&) = &::diy::save<T>
    )
  {
    x.src_gid = src_gid; 
    x.dst_gid = dst_gid;

    const int dst_rank = _assigner.rank(dst_gid);
    if (dst_rank == _async.rank()) {
      _inflight[dst_gid].push_back(std::make_pair(src_gid, x));
      return;
    } 

    _outflight[dst_rank].push_back(x);
  }

  bool dequeue(
      int dst_gid, 
      std::list<std::pair<int, T> > &l, 
      int limit=-1, 
      void (*load)(diy::BinaryBuffer&, T&) = &::diy::load<T>
    )
  {
    typename std::map<int, std::list<std::pair<int, T> > >::iterator it = _inflight.find(dst_gid);
    if (it == _inflight.end() || it->second.empty()) {
      return false;
    }
    else if (limit < 0 || limit >= it->second.size())  // no limit, or the limit is greater than the size
    {
      l = it->second;
      it->second.clear();
      _inflight.erase(it);
      return true;
    } else {
      for (int i=0; i<limit; i++) {
        l.push_back(it->second.front());
        it->second.pop_front();
      }
      return true;
    }
  }

  bool iexchange() // only called in main thread after sync, <del>so no locks required</del>
  {
    if (_status_exchange.global_all_done()) return false;
    _status_exchange.iexchange();

    // incoming
    diy::MemoryBuffer bb;
    int src_rank;
    while (_async.iprobe(src_rank, 0, bb)) {
      std::vector<T> msgs; 
      diy::load(bb, msgs);
      // fprintf(stderr, "in #msgs=%d, size=%d\n", msgs.size(), bb.size());
      
      for (typename std::vector<T>::iterator it = msgs.begin(); it != msgs.end(); it ++) 
      {
        _inflight[it->dst_gid].push_back(std::make_pair(it->src_gid, *it));
      }
    }

    // outgoing
    if (!_async.ready_to_send()) return true;
    std::list<typename std::map<int, std::vector<T> >::iterator> to_remove; 

    for (typename std::map<int, std::vector<T> >::iterator it = _outflight.begin(); 
      it != _outflight.end(); it ++) 
    {
      std::vector<T> &msgs = it->second;
      if (msgs.size() == 0) continue; 
      
      if (!_async.ready_to_send()) 
        break; // avoid serialization cost if comm is not ready

      diy::MemoryBuffer bb;
      diy::save(bb, msgs);
      // fprintf(stderr, "out #msgs=%d, size=%d\n", msgs.size(), bb.size());

      const int dst_rank = it->first; 
      bool succ = _async.isend(dst_rank, 0, bb);
      if (succ) {
        to_remove.push_back(it);
        msgs.clear();
      } // FIXME: else? 
    }

    for (typename std::list<typename std::map<int, std::vector<T> >::iterator>::iterator it = to_remove.begin();
      it != to_remove.end(); it ++)
      _outflight.erase(*it);

    return true;
  }

  bool iexchange(double& itime) // only called in main thread after sync, <del>so no locks required</del>
  {
    if (_status_exchange.global_all_done()) return false;
    _status_exchange.iexchange();

    double start = MPI_Wtime();
    // incoming
    diy::MemoryBuffer bb;
    int src_rank;
    while (_async.iprobe(src_rank, 0, bb)) {
    //if (_async.iprobe(src_rank, 0, bb)) {
      std::vector<T> msgs; 
      diy::load(bb, msgs);
      // fprintf(stderr, "in #msgs=%d, size=%d\n", msgs.size(), bb.size());
      
      for (typename std::vector<T>::iterator it = msgs.begin(); it != msgs.end(); it ++) 
      {
        _inflight[it->dst_gid].push_back(std::make_pair(it->src_gid, *it));
      }
    }

    // outgoing
    if (!_async.ready_to_send()) return true;
    std::list<typename std::map<int, std::vector<T> >::iterator> to_remove; 

    for (typename std::map<int, std::vector<T> >::iterator it = _outflight.begin(); 
      it != _outflight.end(); it ++) 
    {
      std::vector<T> &msgs = it->second;
      if (msgs.size() == 0) continue; 
      
      if (!_async.ready_to_send()) 
        break; // avoid serialization cost if comm is not ready

      diy::MemoryBuffer bb;
      diy::save(bb, msgs);
      // fprintf(stderr, "out #msgs=%d, size=%d\n", msgs.size(), bb.size());

      const int dst_rank = it->first; 
      bool succ = _async.isend(dst_rank, 0, bb);
      if (succ) {
        to_remove.push_back(it);
        msgs.clear();
      } // FIXME: else? 
    }

    for (typename std::list<typename std::map<int, std::vector<T> >::iterator>::iterator it = to_remove.begin();
      it != to_remove.end(); it ++)
      _outflight.erase(*it);

    double end = MPI_Wtime();
    itime += end-start;

    return true;
  }

protected:
  std::map<int, std::vector<T> > _outflight; // key=dst_rank, val=vec(AsyncMsg)
  std::map<int, std::list<std::pair<int, T> > > _inflight; // key=dst_gid, val=list(<src_gid, userMsg>)

  AsyncComm _async;
  const diy::Assigner& _assigner;
  StatusExchange _status_exchange;

};

#endif

