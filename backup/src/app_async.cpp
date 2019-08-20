#include <cassert>
#include <boost/foreach.hpp>
#include "app_async.h"
#include "cons_kdtree.h"

#define NUM_LIMIT 30
#define MAX_NUM_POINTS 100000000

struct VerifyBlockAsync {
  VerifyBlockAsync(CPTApp_Async& app_) :
    app(app_) {}

  void operator()(void* b_, const diy::Master::ProxyWithLink& cp, void*) const {
    Block &b = *static_cast<Block*>(b_);
    const int rank = cp.master()->communicator().rank();

    int count = 0;
    for (size_t i = 0; i < b.particles.size(); i ++) {
      bool flag = true;
      for (int j = 0; j < app.num_dims(); j ++) {
        if (b.particles[i].coords[j] < b.glb[j] || b.particles[i].coords[j] > b.gub[j]) {
          flag = false;
          break;
        }
      }

      if (flag) count ++;
    }
 
    if (count != b.particles.size()) {
      fprintf(stderr, "rank = %d, %d (%d - %d) particles outside the block.\n", 
        rank, (int)b.particles.size() - count, (int)b.particles.size(), count);
      assert (false);
    }
  }

  CPTApp_Async& app;
};

struct TraceBlockAsync { // functor for regular particle tracing
  TraceBlockAsync(CPTApp_Async& app_) :
    app(app_) {}

  void operator()(void *b_, const diy::Master::ProxyWithLink &cp, void*) const {
    Block &b = *static_cast<Block*>(b_);
    AsyncMaster<Message>& async_master = app.asyncMaster();
    const int rank = cp.master()->communicator().rank(), gid = cp.gid();

    // initialize particles
    if (b.num_particles_initialized == -1) {
      app.initialize_particles(b, b.particles);
      b.num_particles_initialized = b.particles.size();
      b.num_particles_traced += b.pending_particles.size();
    }
    
    // dequeue incoming messages
    std::list<std::pair<int, Message> > inmsg;
    async_master.dequeue(gid, inmsg);
    for (std::list<std::pair<int, Message> >::iterator it = inmsg.begin(); 
      it != inmsg.end(); it ++) {
      const Message& msg = it->second; 

      // unfinished particles
      if (msg.unfinished_particles.size() > 0)
        b.particles.insert(b.particles.end(), msg.unfinished_particles.begin(), msg.unfinished_particles.end());

      // finished particles
      if (msg.num_finished_particles > 0) {
        //fprintf(stderr, "[%d] blk %d done, init=%d, done=%d,\n", rank, b.gid, 
        //      (int)b.num_particles_initialized, (int)b.num_particles_finished);
        b.num_particles_finished += msg.num_finished_particles; // pure number of finished particles in this block
        if (b.num_particles_finished == b.num_particles_initialized) {
          fprintf(stderr, "[%d] blk %d done, init=%d, done=%d,\n", rank, b.gid, 
              (int)b.num_particles_initialized, (int)b.num_particles_finished);
          async_master.inc_local_blocks_done();
        }
      }
    }
    inmsg.clear();

    // trace particles
    std::map<int, std::vector<Particle> > unfinished_particles, finished_particles;
    const int limit = NUM_LIMIT; // limit particles for more frequent exchange
    const int count = std::min(limit, (int)b.particles.size());
    std::vector<Particle> working_particles;
    working_particles.insert(working_particles.end(), b.particles.begin(), b.particles.begin() + count);
    b.particles.erase(b.particles.begin(), b.particles.begin() + count);
    if (count > 0) {
      app.trace_particles(b, working_particles, unfinished_particles, finished_particles);
    }

    // enqueue outgoing messages (unfinished particles) 
    std::map<int, Message> outmsg;
    for (std::map<int, std::vector<Particle> >::iterator it = unfinished_particles.begin(); 
        it != unfinished_particles.end(); it ++) 
      outmsg[it->first].unfinished_particles.swap(it->second);

    for (std::map<int, std::vector<Particle> >::iterator it = finished_particles.begin(); 
        it != finished_particles.end(); it ++)
      outmsg[it->first].num_finished_particles = (it->second).size();
    
    for (std::map<int, Message>::iterator it = outmsg.begin();  it != outmsg.end();  it ++) {
      async_master.enqueue(gid, it->first, it->second);
    }
  }

  CPTApp_Async& app;
};
#if 0
struct TraceBlockAsync_Interval { // functor for particle tracing at intervals
  TraceBlockAsync_Interval(CPTApp_Async& app_) :
    app(app_) {}

  void operator()(void *b_, const diy::Master::ProxyWithLink &cp, void*) const {
    Block &b = *static_cast<Block*>(b_);
    AsyncMaster<Message>& async_master = app.asyncMaster();
    const int rank = cp.master()->communicator().rank(), gid = cp.gid();

    // initialize particles at the very beginning
    if (b.num_particles_initialized == -1) {
      app.initialize_particles(b, b.pending_particles);
      b.num_particles_initialized = b.pending_particles.size();
      b.num_particles_traced += b.pending_particles.size();
    } 
    /*
    // initialize particles after intervals (constrained k-d tree decomposition)
    if (b.after_kdtree) { 
      assert (b.pending_particles.size() == 0);
      b.pending_particles.clear();
      b.pending_particles.insert(b.pending_particles.end(), b.particles.begin(), b.particles.end());
      b.particles.clear();
      b.num_particles_initialized = b.pending_particles.size();
      b.num_particles_finished = 0;
      b.after_kdtree = false;
      async_master.reset();
      if (b.num_particles_initialized == 0)
        async_master.inc_local_blocks_done();
    }
   */ 
    // dequeue incoming messages
    std::list<std::pair<int, Message> > inmsg;
    async_master.dequeue(gid, inmsg);
    for (std::list<std::pair<int, Message> >::iterator it = inmsg.begin(); 
      it != inmsg.end(); it ++) {
      const Message& msg = it->second; 

      // unfinished particles
      if (msg.unfinished_particles.size() > 0)
        b.pending_particles.insert(b.pending_particles.end(), msg.unfinished_particles.begin(), msg.unfinished_particles.end());

      // finished particles
      //if (msg.num_finished_particles > 0) 
      {
        b.num_particles_finished += msg.num_finished_particles; // pure number of finished particles in this block
        //fprintf(stderr, "[%d] blk %d done, init=%d, done=%d,\n", rank, b.gid, 
        //      (int)b.num_particles_initialized, (int)b.num_particles_finished);
        if (b.num_particles_finished == b.num_particles_initialized) { // not the real done number To be commented
          //fprintf(stderr, "[%d] blk %d done, init=%d, done=%d,\n", rank, b.gid, 
          //    (int)b.num_particles_initialized, (int)b.num_particles_finished);
          async_master.inc_local_blocks_done();
        }
      }
    }
    inmsg.clear();

    // trace particles
    std::map<int, std::vector<Particle> > unfinished_particles, finished_particles;
    const int limit = NUM_LIMIT; // limit particles for more frequent exchange
    const int count = std::min(limit, (int)b.pending_particles.size());
    std::vector<Particle> working_particles;
    working_particles.insert(working_particles.end(), b.pending_particles.begin(), b.pending_particles.begin() + count);
    b.pending_particles.erase(b.pending_particles.begin(), b.pending_particles.begin() + count);
    if (count > 0) {
      app.trace_particles_kdtree(b, working_particles, unfinished_particles, finished_particles);
    }

    // enqueue outgoing messages (unfinished particles) 
    std::map<int, Message> outmsg;
    for (std::map<int, std::vector<Particle> >::iterator it = unfinished_particles.begin(); 
        it != unfinished_particles.end(); it ++) 
      outmsg[it->first].unfinished_particles.swap(it->second);

    for (std::map<int, std::vector<Particle> >::iterator it = finished_particles.begin(); 
        it != finished_particles.end(); it ++)
      outmsg[it->first].num_finished_particles = (it->second).size();
    
    for (std::map<int, Message>::iterator it = outmsg.begin();  it != outmsg.end();  it ++) {
      async_master.enqueue(gid, it->first, it->second);
    }
  }

  CPTApp_Async& app;
};
#endif
struct TraceBlockAsync_Interval { // functor for particle tracing at intervals
  TraceBlockAsync_Interval(CPTApp_Async& app_) :
    app(app_) {}

  void operator()(void *b_, const diy::Master::ProxyWithLink &cp, void*) const {
    Block &b = *static_cast<Block*>(b_);
    AsyncMaster<Message>& async_master = app.asyncMaster();
    const int rank = cp.master()->communicator().rank(), gid = cp.gid();

    // initialize particles at the very beginning
    if (b.num_particles_initialized == -1) {
      app.initialize_particles(b, b.pending_particles);
      b.num_particles_initialized = b.pending_particles.size();
      b.num_particles_traced += b.pending_particles.size();
      if (b.num_particles_initialized == 0)
        async_master.inc_local_blocks_done();
      //fprintf(stderr, "init = %d\n", (int)b.num_particles_initialized);
    } 

    // initialize particles after intervals (constrained k-d tree decomposition)
    if (b.after_kdtree) { 
      assert (b.pending_particles.size() == 0);
      //b.pending_particles.clear();
      b.pending_particles.insert(b.pending_particles.end(), b.particles.begin(), b.particles.end());
      b.particles.clear();
      b.num_particles_initialized = b.pending_particles.size();
      b.num_particles_finished = 0;
      b.after_kdtree = false;
      //async_master.reset();
      if (b.num_particles_initialized == 0)
        async_master.inc_local_blocks_done();
    }
    
    // dequeue incoming messages
    std::list<std::pair<int, Message> > inmsg;
    async_master.dequeue(gid, inmsg);
    for (std::list<std::pair<int, Message> >::iterator it = inmsg.begin(); 
      it != inmsg.end(); it ++) {
      const Message& msg = it->second; 

      // unfinished particles
      if (msg.unfinished_particles.size() > 0)
        b.pending_particles.insert(b.pending_particles.end(), msg.unfinished_particles.begin(), msg.unfinished_particles.end());

      // finished particles
      if (msg.num_finished_particles > 0) 
      {
        b.num_particles_finished += msg.num_finished_particles; // pure number of finished particles in this block
        //fprintf(stderr, "[%d] blk %d done, init=%d, done=%d,\n", rank, b.gid, 
        //      (int)b.num_particles_initialized, (int)b.num_particles_finished);
        if (b.num_particles_finished == b.num_particles_initialized) { // not the real done number To be commented
          //fprintf(stderr, "[%d] blk %d done, init=%d, done=%d,\n", rank, b.gid, 
          //    (int)b.num_particles_initialized, (int)b.num_particles_finished);
          async_master.inc_local_blocks_done();
        }
      }
    }
    inmsg.clear();

    // trace particles
    std::map<int, std::vector<Particle> > unfinished_particles;//, finished_particles;
    std::map<int, int> num_finished_particles;
    //for (int i=0; i<app.num_blocks(); i++)
    //  num_finished_particles[i] = 0;

    const int limit = NUM_LIMIT; // limit particles for more frequent exchange
    const int count = std::min(limit, (int)b.pending_particles.size());
    std::vector<Particle> working_particles;
    working_particles.insert(working_particles.end(), b.pending_particles.begin(), b.pending_particles.begin() + count);
    b.pending_particles.erase(b.pending_particles.begin(), b.pending_particles.begin() + count);
    if (count > 0) {
      app.trace_particles_kdtree(b, working_particles, unfinished_particles, num_finished_particles);
    }

    // enqueue outgoing messages (unfinished particles) 
    std::map<int, Message> outmsg;
    for (std::map<int, std::vector<Particle> >::iterator it = unfinished_particles.begin(); 
        it != unfinished_particles.end(); it ++) 
      outmsg[it->first].unfinished_particles.swap(it->second);

    for (std::map<int, int>::iterator it = num_finished_particles.begin(); 
        it != num_finished_particles.end(); it ++)
      outmsg[it->first].num_finished_particles = it->second;
    
    for (std::map<int, Message>::iterator it = outmsg.begin();  it != outmsg.end();  it ++) {
      async_master.enqueue(gid, it->first, it->second);
    }
  }

  CPTApp_Async& app;
};

CPTApp_Async::CPTApp_Async() : 
  CPTApp(),
  _time_work(0), 
  _time_kdtree(0),
  _time_exchange(0),
  _time_trace(0),
  _time_iexchange(0),
  _total_steps(0),
  _round_steps(0),
  _max_workload(0),
  _num_particles(0),
  _num_local_initialized(0),
  _num_local_finished(0)
{
}

CPTApp_Async::~CPTApp_Async()
{
  delete _async_master;
}

void CPTApp_Async::init(int argc, char **argv)
{
  CPTApp::init(argc, argv);
  _async_master = new AsyncMaster<Message>(comm_world(), *_assigner, gids().size());
}

void CPTApp_Async::exec()
{
  if (comm_world_rank() == 0) 
    fprintf(stderr, "running...\n");

  if (is_kd_tree()) {  // parallel particle tracing with constrained k-d tree decomposition
    int n_rounds = 0;
    while (true) {
      double t0 = MPI_Wtime();
      _timestamps.push_back(t0); // wait
      _timecategories.push_back(2);
      _master->foreach(TraceBlockAsync_Interval(*this));
      double t1 = MPI_Wtime();
      _timestamps.push_back(t1); // trace
      _timecategories.push_back(1);
      
      double itime = 0;
      bool iexchange = _async_master->iexchange(itime);
      //double t2 = MPI_Wtime();
      _time_trace += t1 - t0;
      _time_iexchange += itime;

      if (!iexchange) {
        _async_master->reset_state();
        if (n_rounds >= NUM_ROUNDS) break;
        else {
          int num_particles_before = 0;
          BOOST_FOREACH(int gid, gids()) {
            Block &b = block(gid);
            b.num_particles_kd = b.particles.size(); // important
            num_particles_before += (int)b.particles.size();
          }

          //Block &b = block(gids()[0]);
          // make statistic before constrained k-d tree decomposition
          //int num_particles_before = (int)b.particles.size(); // TODO int -> size_t
          std::vector<int> point_num1(comm_world_size());
          MPI_Allgather(&num_particles_before, 1, MPI_INT, point_num1.data(), 1, MPI_INT, comm_world());
          int sum_particles = sum_num(point_num1);

          std::vector<int> wl_num(comm_world_size());
          MPI_Allgather(&_round_steps, 1, MPI_INT, wl_num.data(), 1, MPI_INT, comm_world());
          if (comm_world_rank() == 0) {
            _max_workload += max_num(wl_num);
            fprintf(stderr, "--------------- round %d --------------\n", n_rounds);
            fprintf(stderr, "WORKLOAD:\t max = %d, min = %d, avg = %f, max/avg = %f, var = %f\n", 
              max_num(wl_num), min_num(wl_num), avg_num(wl_num), (float)max_num(wl_num)/avg_num(wl_num), var_num(wl_num));
            if (sum_particles != 0)
              fprintf(stderr, "POINTS Before:\t max = %d, min = %d, sum = %d, avg = %f, max/avg = %f, var = %f\n", 
                max_num(point_num1), min_num(point_num1), sum_num(point_num1), avg_num(point_num1), (float)max_num(point_num1)/avg_num(point_num1), var_num(point_num1));
          }
          if (sum_particles == 0) break;
 
          // constrained k-d tree decomposition for particle redistribution
          MPI_Barrier(comm_world());
          double start = MPI_Wtime();
          _timestamps.push_back(start); // wait
          _timecategories.push_back(2);
          // perform constrained k-d tree decomposition (dynamic particle redistribution)
          float perc = (float)sum_particles / (float)MAX_NUM_POINTS;
          if (perc <= 1) { 
            //fprintf(stderr, "rank = %d, sum_particles = %d, perc = %f\n", comm_world_rank(), sum_particles, perc);
            /*BOOST_FOREACH(int gid, gids()) {
              Block &b = block(gid);
              b.particles.insert(b.particles.end(), b.kd_particles.begin(), b.kd_particles.end());
              b.kd_particles.clear();
            }*/
            pt_cons_kdtree_exchange(*_master, *_assigner, _divisions, space_only() ? 3 : _num_dims, space_only(), _block_size, _ghost_size, _constrained, false, true);
          }
          else { // downsampling particles to improve k-d tree decomposition performance
            assert (false);
            /*
            int count = std::ceil(perc);//, each = num_particles_before / count;
            for (int i = 0; i < count; i ++) {
              BOOST_FOREACH(int gid, gids()) {
                Block &b = block(gid);
                int each = (int)b.num_particles_kd / count;
                if (b.kd_particles.size() != 0) {
                  b.particles.insert(b.particles.end(), 
                    b.kd_particles.begin(), i == (count-1) ? b.kd_particles.end() : b.kd_particles.begin() + each);
                  b.kd_particles.erase(b.kd_particles.begin(), i == (count-1) ? b.kd_particles.end() : b.kd_particles.begin() + each);
                }
              }

              pt_cons_kdtree_exchange(*_master, *_assigner, _divisions, space_only() ? 3 : _num_dims, space_only(), _block_size, _ghost_size, _constrained, false, true);

              BOOST_FOREACH(int gid, gids()) {
                Block &b = block(gid);
                if (b.particles.size() != 0) {
                  b.kdtree_particles.insert(b.kdtree_particles.end(), b.particles.begin(), b.particles.end());
                  b.particles.clear();
                }
              }
            }

            BOOST_FOREACH(int gid, gids()) {
              Block &b = block(gid);
              if (b.kdtree_particles.size() != 0) {
                b.particles.insert(b.particles.end(), b.kdtree_particles.begin(), b.kdtree_particles.end());
                b.kdtree_particles.clear();
              }
            }
            */
          }

          //_master->foreach(VerifyBlockAsync(*this));

          MPI_Barrier(comm_world());
          double end = MPI_Wtime();
          _timestamps.push_back(end); // kdtree
          _timecategories.push_back(3);
          _time_kdtree += end - start;

          // make statistic after constrained k-d tree decomposition
          int num_particles_after = 0;
          BOOST_FOREACH(int gid, gids()) {
            Block &b = block(gid);
            num_particles_after += (int)b.particles.size();
          }
          std::vector<int> point_num2(comm_world_size());
          MPI_Allgather(&num_particles_after, 1, MPI_INT, point_num2.data(), 1, MPI_INT, comm_world());
          if (comm_world_rank() == 0) {
            fprintf(stderr, "POINTS After:\t max = %d, min = %d, avg = %f, max/avg = %f, var = %f\n", 
              max_num(point_num2), min_num(point_num2), avg_num(point_num2), (float)max_num(point_num2)/avg_num(point_num2), var_num(point_num2));
          }
/*
          start = MPI_Wtime();
          _timestamps.push_back(start); // wait
          // restart to trace particles
          BOOST_FOREACH(int gid, gids()) {
            Block &b = block(gid);
            if (b.pending_particles.size() != 0) {
              fprintf(stderr, "b.pending_particles.size() = %d\n", (int)b.pending_particles.size());
              assert (false);
            }
            b.pending_particles.insert(b.pending_particles.end(), b.particles.begin(), b.particles.end());
            b.particles.clear();
            b.num_particles_initialized = b.pending_particles.size();
            b.num_particles_finished = 0;
            //b.after_kdtree = false;
            //_async_master->reset();
            if (b.num_particles_initialized == 0)
              _async_master->inc_local_blocks_done();
          }
          end = MPI_Wtime();
          _timestamps.push_back(end); // trace
          _time_trace += end-start;
*/
        }
        _round_steps = 0;
        n_rounds ++;
      }
    }
  }
  else // asynchronous data-parallel particle tracing without particle redistribution 
  { 
#if 0 
    while (_async_master->iexchange()) {
      double t0 = MPI_Wtime();
      _timestamps.push_back(t0); // wait
      _timecategories.push_back(2);
      _master->foreach(TraceBlockAsync(*this));
      double t1 = MPI_Wtime();
      _timestamps.push_back(t1); // trace
      _timecategories.push_back(1);
    }
#else
    int n_rounds = 0;
    while (true) {
      double t0 = MPI_Wtime();
      _timestamps.push_back(t0); // wait
      _timecategories.push_back(2);
      _master->foreach(TraceBlockAsync_Interval(*this));
      double t1 = MPI_Wtime();
      _timestamps.push_back(t1); // trace
      _timecategories.push_back(1);

      double itime = 0;
      bool iexchange = _async_master->iexchange(itime);
      //double t2 = MPI_Wtime();
      _time_trace += t1 - t0;
      _time_iexchange += itime;

      //_master->foreach(TraceBlockAsync_Interval(*this));
      if (!iexchange) {
        _async_master->reset_state();
        if (n_rounds >= NUM_ROUNDS) break;
        else {
          //Block &b = block(gids()[0]);
          int num_particles_before = 0;
          BOOST_FOREACH(int gid, gids()) {
            Block &b = block(gid);
            num_particles_before += (int)b.particles.size();
          }
          std::vector<int> point_num1(comm_world_size());
          MPI_Allgather(&num_particles_before, 1, MPI_INT, point_num1.data(), 1, MPI_INT, comm_world());
          int sum_particles = sum_num(point_num1);

          // make statistic at the intervals
          std::vector<int> wl_num(comm_world_size());
          MPI_Allgather(&_round_steps, 1, MPI_INT, wl_num.data(), 1, MPI_INT, comm_world());
          if (comm_world_rank() == 0) {
            _max_workload += max_num(wl_num);
            fprintf(stderr, "--------------- round %d --------------\n", n_rounds);
            fprintf(stderr, "WORKLOAD:\t max = %d, min = %d, avg = %f, max/avg = %f, var = %f\n", 
              max_num(wl_num), min_num(wl_num), avg_num(wl_num), (float)max_num(wl_num)/avg_num(wl_num), var_num(wl_num));
          }

          if (sum_particles == 0) break;

          BOOST_FOREACH(int gid, gids()) {
            Block &b = block(gid);
            b.after_kdtree = true;
          }
/*
          double start = MPI_Wtime();
          BOOST_FOREACH(int gid, gids()) {
            Block &b = block(gid);
            if (b.pending_particles.size() != 0) {
              fprintf(stderr, "b.pending_particles.size() = %d\n", (int)b.pending_particles.size());
              assert (false);
            }
            b.pending_particles.insert(b.pending_particles.end(), b.particles.begin(), b.particles.end());
            b.particles.clear();
            b.num_particles_initialized = b.pending_particles.size();
            b.num_particles_finished = 0;
            //b.after_kdtree = false;
            //_async_master->reset();
            if (b.num_particles_initialized == 0)
              _async_master->inc_local_blocks_done();
          }
          double end = MPI_Wtime();
          _time_trace += end-start;
*/
        }
        _round_steps = 0;
        n_rounds ++;
      }
    }
#endif
  }

  // finish particles
  if (comm_world_rank() == 0)
    fprintf(stderr, "finishing...\n");
  // TODO add finish particles
  if (comm_world_rank() == 0)
    fprintf(stderr, "done.\n");
}

void CPTApp_Async::stat()
{
  double  _total_time_work, _total_time_exchange;

  MPI_Reduce(
    &_time_work, 
    &_total_time_work, 
    1, MPI_DOUBLE, MPI_SUM, 0, comm_world()
  ); 

  MPI_Reduce(
    &_time_exchange,
    &_total_time_exchange,
    1, MPI_DOUBLE, MPI_SUM, 0, comm_world()
  );

  if (_num_particles == 0) {
    BOOST_FOREACH(int gid, gids()) {
      Block &b = block(gid);
      _num_particles += (int)b.num_particles_traced;
    }
  }

  double time_purerun = _timestamp_finish - _timestamp_init - _time_kdtree,
         time_all     = _timestamp_finish - _timestamp_init + _time_io;

  double time_itrace = _time_trace + _time_iexchange, total_time_itrace;
  double time_iidle = time_purerun - time_itrace, total_time_iidle;

  double time_trace = _time_trace, total_time_trace; // not in use
  double time_idle = time_purerun - time_trace, total_time_idle; // not in use

  MPI_Reduce(
    &time_itrace,
    &total_time_itrace,
    1, MPI_DOUBLE, MPI_SUM, 0, comm_world()
  );

  MPI_Reduce(
    &time_iidle,
    &total_time_iidle,
    1, MPI_DOUBLE, MPI_SUM, 0, comm_world()
  );

  MPI_Reduce( //
    &time_trace,
    &total_time_trace,
    1, MPI_DOUBLE, MPI_SUM, 0, comm_world()
  );

  MPI_Reduce( // 
    &time_idle,
    &total_time_idle,
    1, MPI_DOUBLE, MPI_SUM, 0, comm_world()
  );

  double avg_time_itrace = total_time_itrace/num_blocks(),
         avg_time_iidle = total_time_iidle/num_blocks();

  int _total_num_particles;
  MPI_Reduce(&_num_particles, &_total_num_particles, 1, MPI_INT, MPI_SUM, 0, comm_world());

  std::vector<int> each_num(comm_world_size());
  MPI_Gather(&_total_steps, 1, MPI_INT, each_num.data(), 1, MPI_INT, 0, comm_world());

  if (comm_world_rank() == 0) {
    //for (int i = 0; i < comm_world_size(); i ++)
    //  fprintf(stderr, "proc #%d: integrated %d steps\n", i, each_num[i]);
    fprintf(stderr, "--------------------------------------\n"); 
    fprintf(stderr, "--------------- SUMMARY --------------\n"); 
    fprintf(stderr, "total_num_particles=\t%d\n", _total_num_particles);
    fprintf(stderr, "max_summary_workload=\t%ld\n", _max_workload); 
    fprintf(stderr, "load_balance_indicat=\t%f\n", (float)_max_workload / (float)avg_num(each_num));
    fprintf(stderr, "max_integrated_steps=\t%d\n", max_num(each_num));
    fprintf(stderr, "min_integrated_steps=\t%d\n", min_num(each_num));
    fprintf(stderr, "avg_integrated_steps=\t%f\n", avg_num(each_num));
    fprintf(stderr, "max/avg_integr_steps=\t%.3f\n", (float)max_num(each_num) / avg_num(each_num));
    fprintf(stderr, "var_integrated_steps=\t%f\n", var_num(each_num));
    fprintf(stderr, "total_time_work=\t%f\n", _total_time_work); 
    fprintf(stderr, "total_time_exchange=\t%f\n", _total_time_exchange);
    fprintf(stderr, "--------------------------------------\n");
    fprintf(stderr, "time_init=\t\t%.3f\n", _timestamp_init - _timestamp_start); 
    fprintf(stderr, "time_run=\t\t%.3f\n", _timestamp_finish - _timestamp_init); 
    fprintf(stderr, "time_total=\t\t%.3f\n", _timestamp_finish - _timestamp_start); 
    fprintf(stderr, "--------------------------------------\n");
    fprintf(stderr, "time_avgtrace=\t\t%.3f\n", total_time_trace/num_blocks());
    fprintf(stderr, "time_avgidle=\t\t%.3f\n", total_time_idle/num_blocks());
    fprintf(stderr, "--------------------------------------\n");
    fprintf(stderr, "time_purerun=\t\t%.3f\t\t%.3f\n", time_purerun, time_purerun/time_all);
    fprintf(stderr, "time_avgitrace=\t\t%.3f\t\t%.3f\n", avg_time_itrace, avg_time_itrace/time_all);
    fprintf(stderr, "time_avgiidle=\t\t%.3f\t\t%.3f\n", avg_time_iidle, avg_time_iidle/time_all);
    fprintf(stderr, "time_io=\t\t%.3f\t\t%.3f\n", _time_io, _time_io/time_all);
    fprintf(stderr, "time_kdtree=\t\t%.3f\t\t%.3f\n", _time_kdtree, _time_kdtree/time_all);
    fprintf(stderr, "time_all=\t\t%.3f\n", time_all);
    fprintf(stderr, "--------------------------------------\n"); 
  }
}

void CPTApp_Async::add_workload()
{
  _round_steps ++;
  _total_steps ++;
}
