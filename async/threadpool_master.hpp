#ifndef _THREADPOOL_MASTER_H
#define _THREADPOOL_MASTER_H

#include <diy/assigner.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/function.hpp>
#include <boost/foreach.hpp>
#include <boost/thread/tss.hpp>
#include <boost/atomic.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/container/vector.hpp> // container of non-copyable objs
#include <vector>
#include <queue>
#include <sched.h>
#include <sys/syscall.h>
#include "status.h"
#include "comm.h"

#if WITH_CXX11
#include "concurrentqueue/concurrentqueue.h"
#include "concurrentqueue/blockingconcurrentqueue.h"
#endif

#if 0
typedef AsyncComm_MY AsyncComm;
#else
typedef AsyncComm_HG AsyncComm;
#endif

template<class BlkMsg, class ProcStat> // T=msg. 
class ThreadPoolMaster {
#if WITH_CXX11
  template <class T>
  using concurrent_queue = moodycamel::ConcurrentQueue<T>;
#endif

public:
  ThreadPoolMaster(MPI_Comm comm, const diy::Assigner& assigner, int num_local_blocks) :
    _comm(comm),
    _assigner(assigner),
    _status_exchange(comm, num_local_blocks),
    _async(comm, 10), 
    _qwork(65536)
  {
    MPI_Comm_size(comm, &_np);
    MPI_Comm_rank(comm, &_rank);
#if WITH_CXX11
    _qsends.resize(_np); 
#else
    for (int i=0; i<_np; i++)
      _qsends.push_back(new boost::lockfree::queue<BlkMsg*>(65536));
#endif
  }
  ~ThreadPoolMaster() {} // FIXME: delete _qsends
  
  uint64_t num_bytes_sent() const {return _async.num_bytes_sent();}

  void inc_local_blocks_done() {_status_exchange.inc_local_blocks_done();}
  void set_local_blocks_done(int n) {_status_exchange.set_local_blocks_done(n);}

  void set_msg_cb(boost::function<void(int, BlkMsg&)> f) {_f = f;} // handle msg (usually in worker thread, heavy work)
  //void set_cycle_cb(boost::function<void()> h) {_h = h;} // called every cycle in the comm thread

  void set_affinity(int processorID) {
#if WITH_AFFINITY
    cpu_set_t cpu_set;
    CPU_ZERO(&cpu_set);
    CPU_SET(processorID, &cpu_set); // the 0th CPU is for comm
   
    pthread_t thread = pthread_self();
    // pid_t pid = syscall(SYS_gettid);
    // int rtn = sched_setaffinity(pid, sizeof(cpu_set_t), &cpu_set);
    int rtn = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpu_set);
    // int cpu = sched_getcpu();
    // fprintf(stderr, "rank=%d, processorID=%d, rtn=%d\n", 
    //     _rank, processorID, rtn);
#endif
  }

  void run(int nthreads)
  {
    boost::thread_group worker_threads;
    for (int tid = 1; tid < nthreads; tid ++)
      worker_threads.create_thread(boost::bind(&ThreadPoolMaster::work, this, tid));

    communicate();
    worker_threads.join_all(); 
  }

  void work(int tid) { // consume (trace particles) and produce (enqueue particles)
    set_affinity(tid);

#if WITH_CXX11
    BlkMsg msg;
    while (!_status_exchange.global_all_done()) {
      if (_qwork.try_dequeue(msg))
        _f(msg.dst_gid, msg);
    }
#else
    while (!_status_exchange.global_all_done()) {
      BlkMsg *m = NULL;
      if (_qwork.pop(m))
        _f(m->dst_gid, *m);
      if (m != NULL) delete m;
    }
#endif
  }

  void communicate() { // produce
    set_affinity(0);
    while (!_status_exchange.global_all_done()) {
      for (int dst_rank = 0; dst_rank < _np; dst_rank ++) {
        if (dst_rank == _rank) continue;
        if (!_async.ready_to_send()) continue;

#if WITH_CXX11
        static const int bulk_size = 2048;
        static BlkMsg msgs_[bulk_size];
        size_t count = _qsends[dst_rank].try_dequeue_bulk(msgs_, bulk_size);
        if (count == 0) continue;
        // fprintf(stderr, "bulk_send=%d\n", count);

        //ProcStat ps; 
        //_pf(ps);
        std::vector<BlkMsg> msgs(msgs_, msgs_ + count);

        diy::MemoryBuffer bb;
        //diy::save(bb, ps);
        diy::save(bb, msgs);
        
        if (!_async.isend(dst_rank, 0, bb)) {
          fprintf(stderr, "rank=%d, Oops\n", _rank);
          assert(false);
          for (int i=0; i<msgs.size(); i++)
            _qsends[dst_rank].enqueue(msgs[i]);
        }
#else
        std::vector<BlkMsg> msgs;
        BlkMsg *m = NULL; 
        while (_qsends[dst_rank]->pop(m)) {
          msgs.push_back(BlkMsg());
          msgs.back().swap(*m);
          delete m;
        }

        if (msgs.size() > 0) {
          // fprintf(stderr, "%d->%d, %d\n", _rank, dst_rank, msgs.size());
          //ProcStat ps; 
          //_pf(ps);

          diy::MemoryBuffer bb;
          //diy::save(bb, ps);
          diy::save(bb, msgs);
          
          if (!_async.isend(dst_rank, 0, bb)) {
            fprintf(stderr, "rank=%d, Oops\n", _rank);
            assert(false);
          }
        }
#endif
      }
      
      // iprobe and irecv
      int count = 0;
      int src_rank; 
      diy::MemoryBuffer bb;
      while (_async.iprobe(src_rank, 0, bb)) {
        // ProcStat ps;  // pstat, the "header"
        // diy::load(bb, ps);
        // _pg(src_rank, ps);
        std::vector<BlkMsg> msgs;
        diy::load(bb, msgs);
        // fprintf(stderr, "%d<-%d, %d\n", _rank, src_rank, msgs.size());
        for (int i=0; i<msgs.size(); i++) {
          enqueue(msgs[i]);
        }
        count ++;
      }

      //if (!_h.empty()) _h(); // called every cycle
      _status_exchange.iexchange();
    }
  }

  void enqueue(int src_gid, int dst_gid, BlkMsg& x) { // called from worker threads, which are both producers and consumers
    x.src_gid = src_gid; 
    x.dst_gid = dst_gid;
    enqueue(x);
  }

protected:
  void enqueue(BlkMsg& msg) {
    const int dst_rank = _assigner.rank(msg.dst_gid);
#if WITH_CXX11
    if (dst_rank == _rank) 
      _qwork.enqueue(msg);
    else 
      _qsends[dst_rank].enqueue(msg);
#else
    BlkMsg *m = new BlkMsg();
    m->swap(msg);
    if (dst_rank == _rank)
      _qwork.push(m);
    else
      _qsends[dst_rank]->push(m);
#endif
  }

protected:
#if WITH_CXX11
  concurrent_queue<BlkMsg> _qwork;
  std::vector<concurrent_queue<BlkMsg> > _qsends; // dst_rank, qsend
#else
  boost::lockfree::queue<BlkMsg*> _qwork;
  std::vector<boost::lockfree::queue<BlkMsg*>*> _qsends;
#endif

  boost::function<void(int, BlkMsg&)> _f;
  //boost::function<bool(int, BlkMsg&)> _g;
  //boost::function<void()> _h;

  MPI_Comm _comm;
  AsyncComm _async;
  StatusExchange _status_exchange;

  const diy::Assigner& _assigner;
  int _np, _rank;
};

#endif
