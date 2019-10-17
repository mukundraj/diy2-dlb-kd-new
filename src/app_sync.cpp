#include <cassert>
#include <boost/foreach.hpp>
#include "common/utils.h"
#include "app_sync.h"
#include "cons_kdtree.h"

#define DEBUG 0
#define STORE_PARTICLES 0

struct TraceBlockSync { // functor for regular particle tracing
  TraceBlockSync(CPTApp_Sync& app_) :
    app(app_) {}

  void operator()(void *b_, const diy::Master::ProxyWithLink &cp, void*) const {
    Block &b = *static_cast<Block*>(b_);
    const int rank = cp.master()->communicator().rank(), gid = cp.gid(); 
    std::map<int, std::vector<Particle> > unfinished_particles, finished_particles;
    std::vector<Particle> incoming_particles;

    // initialize seed particles first time
    if (b.num_particles_initialized == -1) {
      app.initialize_particles(b, incoming_particles);
      b.num_particles_initialized = incoming_particles.size();
      app._local_init += b.num_particles_initialized;
    }
    
    // dequeue vectors of endpoints, add to seed particles
    std::vector<int> in;
    cp.incoming(in);
    for (int i = 0; i < in.size(); i ++) {
      if (cp.incoming(in[i]).buffer.size() > 0) {
        std::vector<Particle> ps;
        cp.dequeue(in[i], ps);
        incoming_particles.insert(incoming_particles.end(), ps.begin(), ps.end());
      }
    }

    // add unfinished particles in the previous round
    for (int i = 0; i < b.particles.size(); i ++)
      incoming_particles.push_back(b.particles[i]);
    b.particles.clear();

    // trace particles
    if (incoming_particles.size() > 0) {
      app.trace_particles(b, incoming_particles, unfinished_particles, finished_particles);
      incoming_particles.clear();
    
      for (std::map<int, std::vector<Particle> >::iterator it = finished_particles.begin(); it != finished_particles.end(); it ++) {
        b.num_particles_finished += it->second.size();
      }
    }

    // enqueue the vectors of endpoints
    for (std::map<int, std::vector<Particle> >::iterator it = unfinished_particles.begin(); it != unfinished_particles.end(); it ++) {
      diy::BlockID bid; 
      bid.gid = it->first; 
      bid.proc = app.rank(bid.gid);

      if (bid.gid == b.gid) { 
        b.particles.insert(b.particles.end(), (it->second).begin(), (it->second).end());
      }
      else {
        bool flag = false;
        for (int i=0; i<cp.link()->size(); i++) {
          if (bid.gid == cp.link()->target(i).gid) {
            flag = true;
            break;
          }
        }
        if (!flag) {
          fprintf(stderr, "ERROR: bid is not in cp.link neighbors\n");
          assert (false);
        }
        cp.enqueue(bid, it->second);
      }
    }

    // stage all_reduce of total initialized and total finished particle traces
    cp.all_reduce((int)b.num_particles_initialized, std::plus<int>());
    cp.all_reduce((int)b.num_particles_finished, std::plus<int>());
  }

  CPTApp_Sync& app;
};

size_t _local_unfinished_particles, _local_num_particles;

struct TraceBlockRound { // functor for doing a round of particle tracing in epoch pred load bal
  TraceBlockRound(CPTApp_Sync& app_) :
    app(app_) {}

  void operator()(void *b_, const diy::Master::ProxyWithLink &cp) const {
    Block &b = *static_cast<Block*>(b_);
    const int rank = cp.master()->communicator().rank(), gid = cp.gid(); 
    std::map<int, std::vector<Particle> > unfinished_particles, finished_particles;
    std::vector<Particle> incoming_particles;

    // diy::Link* l = cp.link();
    // fprintf(stderr, "Link size : %d\n", l->size());
    // for (int i = 0; i < l->size(); ++i) {
    //   fprintf(stderr, "lt %d:-%d,", b.gid, l->target(i));
    // }
    // fprintf(stderr, "\n");
   
    // dequeue vectors of endpoints, add to seed particles
    std::vector<int> in;
    cp.incoming(in);
    for (int i = 0; i < in.size(); i ++) {
      if (cp.incoming(in[i]).buffer.size() > 0) {
        std::vector<Particle> ps;
        cp.dequeue(in[i], ps);
        incoming_particles.insert(incoming_particles.end(), ps.begin(), ps.end());
        for (size_t ii=0; ii<ps.size(); ii++){
          fprintf (stderr, "[%d (%d)] ", b.gid, ps[ii].id);
        }
        fprintf (stderr, "\n ");
      }

    }

    _local_num_particles = incoming_particles.size();

    // add unfinished particles in the previous round
    for (int i = 0; i < b.particles.size(); i ++)
      incoming_particles.push_back(b.particles[i]);
    b.particles.clear();


    // trace particles
    if (incoming_particles.size() > 0) {
      app.trace_particles_core(b, incoming_particles, unfinished_particles, finished_particles, cp);
      incoming_particles.clear();
    
      for (std::map<int, std::vector<Particle> >::iterator it = finished_particles.begin(); it != finished_particles.end(); it ++) {
        b.num_particles_finished += it->second.size();
      }
    }


    _local_unfinished_particles = 0;
    // enqueue the vectors of endpoints
    for (std::map<int, std::vector<Particle> >::iterator it = unfinished_particles.begin(); it != unfinished_particles.end(); it ++) {
      diy::BlockID bid; 
      bid.gid = it->first; 
      bid.proc = app.rank(bid.gid);

      _local_unfinished_particles += it->second.size();

      if (bid.gid == b.gid) { 
        b.particles.insert(b.particles.end(), (it->second).begin(), (it->second).end());
      }
      else {
        bool flag = false;
        for (int i=0; i<cp.link()->size(); i++) {
          if (bid.gid == cp.link()->target(i).gid) {
            flag = true;
            break;
          }


        }
        if (!flag) {
          fprintf(stderr, "ERROR: bid is not in cp.link neighbors\n");
          assert (false);
        }
        for (size_t ii=0; ii<it->second.size(); ii++){
                  fprintf (stderr, "{%d->%d (%d)} ", b.gid, bid, it->second[ii].id);
        }
         fprintf (stderr, "\n");
        
        cp.enqueue(bid, it->second);
        

      }
    }

    // stage all_reduce of total initialized and total finished particle traces
    cp.all_reduce((int)b.num_particles_initialized, std::plus<int>());
    cp.all_reduce((int)b.num_particles_finished, std::plus<int>());
  }

  CPTApp_Sync& app;
  
};


struct VerifyBlockSync {
  VerifyBlockSync(CPTApp_Sync& app_) :
    app(app_) {}

  void operator()(void* b_, const diy::Master::ProxyWithLink& cp, void*) const {
    Block &b = *static_cast<Block*>(b_);
    const int rank = cp.master()->communicator().rank();

    int count = 0;
    for (size_t i = 0; i < b.particles.size(); i ++) {
      bool flag = true;
      for (int j = 0; j < app.num_dims(); j ++) {
        if (b.particles[i].coords[j] < b.glb[j] || b.particles[i].coords[j] >= b.gub[j]) {
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

  CPTApp_Sync& app;
};



CPTApp_Sync::CPTApp_Sync() : 
  CPTApp(),
  _time_work(0), 
  _time_kdtree(0),
  _time_exchange(0),
  _time_trace(0),
  _local_init(0),
  _local_done(0),
  _total_steps(0),
  _round_steps(0),
  _epoch_steps(0), 
  _max_workload(0),
  _num_particles(0),
  _pred_mismatch(0),
  _total_pred_mismatch(0)
{
}

void CPTApp_Sync::init(int argc, char **argv)
{
  CPTApp::init(argc, argv);
}

void CPTApp_Sync::exec()
{
  if (comm_world_rank() == 0) 
    fprintf(stderr, "running...\n");

  int count = 0;
  if (is_kd_tree()) {  // parallel particle tracing with constrained k-d tree decomposition
    int epoch_ctr = 0;
    while (true) { // each iteration is an Epoch
      // double t0 = MPI_Wtime();
      //_timestamps.push_back(t0); // wait
      //_timecategories.push_back(2);

#if STORE_PARTICLES
      std::vector<Particle> particles_first, particles_second;
#endif

       BOOST_FOREACH(int gid, gids()) {
        Block &b = block(gid);

                // initialize particles at the first time execution
        if (b.num_particles_initialized == -1) {
          initialize_particles(b, b.particles);
          b.num_particles_initialized = b.particles.size();
          b.num_particles_traced += b.particles.size();
          _local_init += b.num_particles_initialized;
          #if STORE_PARTICLES
                    gather_store_particles(count, b.particles); count ++;
                    gather_store_cores(b); // core bounds for split
          #endif
        }




        
      }



      // prediction and generate ghost particles (not for baseline case)


      // compute kd-tree based on currently stored points
      pt_cons_kdtree_exchange(*_master, *_assigner, _divisions, space_only() ? 3 : _num_dims, space_only(), _block_size, _ghost_size, _constrained, false, false);



      // if prediction case, then filter out ghost particles


      _master->foreach([&](Block* b, const diy::Master::ProxyWithLink& cp) { 
      RCLink*  link      = static_cast<RCLink*>(cp.link());
       fprintf(stderr, "gid %d, Lsize %d\n", cp.gid(), link->size());
       for (int i = 0; i < b->nbr_gids.size(); ++i)
        { fprintf(stderr, " %d brgid %d (%f %f, %f %f, %f %f)\n", cp.gid(), b->nbr_gids[i], b->nbr_bounds[i*6+0], b->nbr_bounds[i*6+1], b->nbr_bounds[i*6+2], b->nbr_bounds[i*6+3], b->nbr_bounds[i*6+4], b->nbr_bounds[i*6+5]);
        }
      });


      _local_init_epoch = 0;
       BOOST_FOREACH(int gid, gids()) {
        Block &b = block(gid);
        _local_init_epoch += b.particles.size();


        int gst[4], gsz[4], lst[4], lsz[4];
        b.get_ghost_load_st_sz(num_dims(), gst, gsz, lst, lsz);
        float clb[4], cub[4]; // core_start and core_size
        b.get_core_st_sz(num_dims(), clb, cub);
        fprintf(stderr, "gid: %d, (%f %f) (%f %f) (%f %f)// (%d %d %d) (%d %d %d) // (%d %d %d) (%d %d %d) \n", gid, clb[0], cub[0], clb[1], cub[1], clb[2], cub[2], 
          gst[0], gst[1], gst[2], gsz[0], gsz[1], gsz[2],
          lst[0], lst[1], lst[2], lsz[0], lsz[1], lsz[2]);

      }
      _local_done_epoch = 0;

      int ctr = 0;
      while(true){ // epoch loop: each iteration is a round
          int init_epoch = 0, done_epoch = 0;
          size_t total_particles = 0, total_unfinished_particles;

         

          // _master->foreach(TraceBlockRound(*this));
          _master->foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
                    {
                        // fprintf(stderr, "ctr %d, rank %d, nparticles %ld\n", ctr, comm_world_rank(), b->particles.size() );
                        TraceBlockRound tbr(*this);
                        tbr(b, cp);

          });

          // bool remote = true;
          _master->exchange();
            

           MPI_Barrier(comm_world());
            fprintf(stderr, "ctr %d done\n", ctr);
            MPI_Barrier(comm_world());

           MPI_Allreduce(&_local_num_particles, &total_particles, 1, MPI_LONG, MPI_SUM, comm_world());
           MPI_Allreduce(&_local_unfinished_particles, &total_unfinished_particles, 1, MPI_LONG, MPI_SUM, comm_world());


          // check if epoch is done
          MPI_Allreduce(&_local_init_epoch, &init_epoch, 1, MPI_INT, MPI_SUM, comm_world());
          MPI_Allreduce(&_local_done_epoch, &done_epoch, 1, MPI_INT, MPI_SUM, comm_world());

          if (comm_world_rank() == 0){
            fprintf(stderr, "init_epoch %d, done_epoch %d, _local_num_particles %ld, total_particles (tot incoming at st o rnd) %ld, total_unfinished_particles %ld\n", init_epoch, done_epoch, _local_num_particles, total_particles, total_unfinished_particles);
          }
         
          ctr++;
          if (ctr==10 ) break;



          std::vector<int> round_wl_num(comm_world_size());
          MPI_Allgather(&_round_steps, 1, MPI_INT, round_wl_num.data(), 1, MPI_INT, comm_world());
          if (comm_world_rank() == 0)
            _round_balance.push_back((float)max_num(round_wl_num)/avg_num(round_wl_num));

          _round_steps = 0;

          if (init_epoch == done_epoch){ // ONLY FOR TESTING: inner loop break condition (breaks from epoch)
            break; 
          } 
      } // end while : round loop

      epoch_ctr++;
      MPI_Barrier(comm_world());
      if (comm_world_rank() == 0){

        fprintf(stderr, "Epochs completed: %d\n", epoch_ctr);
      }


//     MPI_Barrier(comm_world());
//     double _time_pred_start = MPI_Wtime();   
     
//     BOOST_FOREACH(int gid, gids()) {
//         Block &b = block(gid);

//         if (pred_val()>0){
//           // predict workload
//           trace_particles_kdtree_predict(b, pred_val());
//         }

//       }


//       MPI_Barrier(comm_world());
//       double kd_start = MPI_Wtime();
//       _time_prediction += kd_start - _time_pred_start;
//       int num_particles_before = 0;

//        BOOST_FOREACH(int gid, gids()) {
//         Block &b = block(gid);

//         // load balance 
//         pt_cons_kdtree_exchange(*_master, *_assigner, _divisions, space_only() ? 3 : _num_dims, space_only(), _block_size, _ghost_size, _constrained, false, false); 

//       }

//       MPI_Barrier(comm_world());
//       double end = MPI_Wtime();
//       _timestamps.push_back(end); // kdtree
//       _timecategories.push_back(3);
//       _time_kdtree += end - kd_start;


//     double t0 = MPI_Wtime();
//     BOOST_FOREACH(int gid, gids()) {
//         Block &b = block(gid);

//         // trace particles 
//         std::map<int, std::vector<Particle> > unfinished_particles, finished_particles;
//         std::vector<Particle> working_particles;
//         working_particles.insert(working_particles.end(), b.particles.begin(), b.particles.end());
//         b.particles.clear();
//         if (working_particles.size() > 0) {
// #if STORE_PARTICLES
//           trace_particles_video(b, working_particles, particles_first, particles_second);
// #else
//           trace_particles_kdtree(b, working_particles, unfinished_particles, finished_particles);
// #endif
//         }

//         num_particles_before += (int)b.particles.size();
//       }

//       double t1 = MPI_Wtime();
//       _timestamps.push_back(t1); // trace
//       _timecategories.push_back(1);

//       _time_trace += t1 - t0;

      // check if finished particle tracing
      int init, done; 
      MPI_Allreduce(&_local_init, &init, 1, MPI_INT, MPI_SUM, comm_world());
      MPI_Allreduce(&_local_done, &done, 1, MPI_INT, MPI_SUM, comm_world());

      std::vector<int> wl_num(comm_world_size());
      MPI_Allgather(&_epoch_steps, 1, MPI_INT, wl_num.data(), 1, MPI_INT, comm_world());

      if (comm_world_rank() == 0) 
        _max_workload += max_num(wl_num);
      if (comm_world_rank() == 0)
        _balance.push_back((float)max_num(wl_num)/avg_num(wl_num));

#if DEBUG
      if (comm_world_rank() == 0) {
        //_max_workload += max_num(wl_num);
        fprintf(stderr, "--------------- round %d --------------\n", _nt_loops);
        fprintf(stderr, "INIT=%d, DONE=%d\n", init, done); // round
        fprintf(stderr, "WORKLOAD:\t max = %d, min = %d, avg = %f, max/avg = %f, var = %f\n", 
          max_num(wl_num), min_num(wl_num), avg_num(wl_num), (float)max_num(wl_num)/avg_num(wl_num), var_num(wl_num));
        //_balance.push_back((float)max_num(wl_num)/avg_num(wl_num));
      }
#endif

      _epoch_steps = 0;
      _nt_loops ++;

#if STORE_PARTICLES
      BOOST_FOREACH(int gid, gids()) {
        Block &b = block(gid);
        gather_store_particles(count, particles_first); count ++;
        gather_store_particles(count, particles_second); count ++;
        gather_store_particles(count, b.particles); count ++;
        //gather_store_cores(b); // core bounds for split
      }
#endif

      // if (init == done && done != 0) break; // main loop (outer loop) break condition
      if (init == done) break; // main loop (outer loop) break condition


#if DEBUG
      std::vector<int> point_num1(comm_world_size());
      MPI_Allgather(&num_particles_before, 1, MPI_INT, point_num1.data(), 1, MPI_INT, comm_world());
      int sum_particles = sum_num(point_num1);
      if (comm_world_rank() == 0) {
        if (sum_particles != 0) 
          fprintf(stderr, "POINTS Before:\t max = %d, min = %d, sum = %d, avg = %f, max/avg = %f, var = %f\n", 
              max_num(point_num1), min_num(point_num1), sum_num(point_num1), avg_num(point_num1), (float)max_num(point_num1)/avg_num(point_num1), var_num(point_num1));
      }

      if (sum_particles == 0) break;
#endif
      MPI_Barrier(comm_world());
      double start = MPI_Wtime();
      _timestamps.push_back(start); // wait
      _timecategories.push_back(2);
      // perform constrained k-d tree decomposition (dynamic particle redistribution)
      // pt_cons_kdtree_exchange(*_master, *_assigner, _divisions, space_only() ? 3 : _num_dims, space_only(), _block_size, _ghost_size, _constrained, false, false); 
      //_master->foreach(VerifyBlockSync(*this));

      // MPI_Barrier(comm_world());
      // double end = MPI_Wtime();
      // _timestamps.push_back(end); // kdtree
      // _timecategories.push_back(3);
      // _time_kdtree += end - start;
#if DEBUG
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
#endif

#if STORE_PARTICLES
      BOOST_FOREACH(int gid, gids()) {
        Block &b = block(gid);
        gather_store_particles(count, b.particles); count ++;
        gather_store_cores(b); // core bounds for split
      }
#endif
    } // end while : epoch loop
  }
  else {  // data-parallel particle tracing without redistribution (TOD using MPI_Isend and MPI_Irecv)
    while (1) {
      double t0 = MPI_Wtime();
      //_timestamps.push_back(t0); // wait
      //_timecategories.push_back(2);

      // _master->foreach(TraceBlockSync(*this)); // trace particles
      _master->foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
                    {
                        TraceBlockSync(*this);
        });

      double t1 = MPI_Wtime();
      _timestamps.push_back(t1); // trace
      _timecategories.push_back(1);

      _time_trace += t1 - t0;

      MPI_Barrier(comm_world());
      double start = MPI_Wtime();
      _timestamps.push_back(start); // wait
      _timecategories.push_back(2);
      _master->exchange();  
      MPI_Barrier(comm_world());
      double end = MPI_Wtime();
      //_timestamps.push_back(end);
      //_timecategories.push_back(2);
      _timestamps.push_back(end); // exchange (similar to kdtree)
      _timecategories.push_back(3);

      _time_kdtree += end - start;

      int init, done;
      for (int i=0; i<gids().size(); i++) {
        init = _master->proxy(i).get<int>();
        done = _master->proxy(i).get<int>();
      }

      //MPI_Allreduce(&_local_init, &init, 1, MPI_INT, MPI_SUM, comm_world());
      //MPI_Allreduce(&_local_done, &done, 1, MPI_INT, MPI_SUM, comm_world());

      std::vector<int> wl_num(comm_world_size());
      MPI_Allgather(&_round_steps, 1, MPI_INT, wl_num.data(), 1, MPI_INT, comm_world());

      if (comm_world_rank() == 0) 
        _max_workload += max_num(wl_num);
      if (comm_world_rank() == 0)
        _balance.push_back((float)max_num(wl_num)/avg_num(wl_num));

#if DEBUG
      if (comm_world_rank() == 0) {
        //_max_workload += max_num(wl_num);
        fprintf(stderr, "--------------- round %d --------------\n", _nt_loops);
        fprintf(stderr, "INIT=%d, DONE=%d\n", init, done); // round
        fprintf(stderr, "WORKLOAD:\t max = %d, min = %d, avg = %f, max/avg = %f, var = %f\n", 
          max_num(wl_num), min_num(wl_num), avg_num(wl_num), (float)max_num(wl_num)/avg_num(wl_num), var_num(wl_num));
        //_balance.push_back((float)max_num(wl_num)/avg_num(wl_num));
      }
#endif
      _round_steps = 0;
      _nt_loops ++;

      if (init == done && done != 0) // main loop break for second non kdtree condition
        break;

    }
  }

  // finish particles
  if (comm_world_rank() == 0)
    fprintf(stderr, "finishing...\n");
  // TODO add finish particles
  if (comm_world_rank() == 0)
    fprintf(stderr, "done.\n");
}

void CPTApp_Sync::stat()
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

  MPI_Reduce(
    &_pred_mismatch,
    &_total_pred_mismatch,
    1, MPI_INT, MPI_SUM, 0, comm_world()
  );

  if (_num_particles == 0) {
    BOOST_FOREACH(int gid, gids()) {
      Block &b = block(gid);
      _num_particles += (int)b.num_particles_traced;
    }
  }

  int _total_num_particles;
  MPI_Reduce(&_num_particles, &_total_num_particles, 1, MPI_INT, MPI_SUM, 0, comm_world());

  double time_purerun = _timestamp_finish - _timestamp_init - _time_kdtree,
         time_all     = _timestamp_finish - _timestamp_init + _time_io;

  double time_trace = _time_trace, total_time_trace; 
  double time_idle = time_purerun - time_trace, total_time_idle; 

  MPI_Reduce( 
    &time_trace,
    &total_time_trace,
    1, MPI_DOUBLE, MPI_SUM, 0, comm_world()
  );

  MPI_Reduce( 
    &time_idle,
    &total_time_idle,
    1, MPI_DOUBLE, MPI_SUM, 0, comm_world()
  );

  double avg_time_trace = total_time_trace/num_blocks(),
         avg_time_idle = total_time_idle/num_blocks();

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
    fprintf(stderr, "time_avgtrace=\t\t%.3f\t\t%.3f\n", avg_time_trace, avg_time_trace/time_all);
    fprintf(stderr, "time_avgidle=\t\t%.3f\t\t%.3f\n", avg_time_idle, avg_time_idle/time_all);
    fprintf(stderr, "time_io=\t\t%.3f\t\t%.3f\n", _time_io, _time_io/time_all);
    fprintf(stderr, "time_redistrib=\t\t%.3f\t\t%.3f\n", _time_kdtree, _time_kdtree/time_all);
    fprintf(stderr, "time_all=\t\t%.3f\n", time_all);
    fprintf(stderr, "_time_trace=\t\t%.3f\n", _time_trace);
    fprintf(stderr, "_time_kdtree=\t\t%.3f\n", _time_kdtree);
    fprintf(stderr, "_time_prediction=\t\t%.3f\n", _time_prediction);
    fprintf(stderr, "_total_pred_mismatch=\t\t%d\n", _total_pred_mismatch);
    fprintf(stderr, "--------------------------------------\n"); 
  }
} 

void CPTApp_Sync::gather_store_cores(const Block& b)
{
  std::vector<std::string> local_entities(1);

  PBBound bound;
  bound.set_x_min(b.core_bounds.min[0]);
  bound.set_x_max(b.core_bounds.max[0]);
  bound.set_y_min(b.core_bounds.min[1]);
  bound.set_y_max(b.core_bounds.max[1]);

  bound.SerializeToString(&local_entities[0]);

  assert(_num_dims == 2);
  int total_num;
  std::vector<std::string> all_entities;
  gather_protobuf(1, local_entities, total_num, all_entities, 0, comm_world());
        
  assert(total_num == comm_world_size());
  if (comm_world_rank() == 0) {
    fprintf(stderr, "--------------- round %d --------------\n", _nt_loops);
    for (int i = 0; i < total_num; i ++) {
      PBBound bound;
      bound.ParseFromString(all_entities[i]);
      fprintf(stderr, "%d, [%f, %f] -> [%f, %f]\n", i, bound.x_min(),
        bound.y_min(), bound.x_max(), bound.y_max());
    }
  }
}

void CPTApp_Sync::gather_store_particles(const int count, const std::vector<Particle>& particles)
{
#if 0
  int local_num = (int)particles.size();
  std::vector<std::string> local_entities(local_num);
  int n = 0;

  BOOST_FOREACH(Particle p, particles) {
    PBParticle particle;
    particle.set_x(p[0]);
    particle.set_y(p[1]);
    if (_num_dims >= 3) particle.set_z(p[2]);
    if (_num_dims == 4) particle.set_t(p[3]);

    particle.SerializeToString(&local_entities[n++]);
  }

  assert(_num_dims == 2);
  int total_num;
  std::vector<std::string> all_entities;
  gather_protobuf(local_num, local_entities, total_num, all_entities, 0, comm_world());
        
  if (comm_world_rank() == 0) {
    std::stringstream stream;  
    stream << "video/particles/" << count << ".csv"; 

    std::ofstream ofile;  
    ofile.open(stream.str());

    for (int i = 0; i < total_num; i ++) {
      PBParticle particle;
      particle.ParseFromString(all_entities[i]);
      ofile << particle.x() << ":" << particle.y() << ",";
    }
    ofile.close();
  }
#else
  std::stringstream stream;  
  stream << "video/particles/" << count << "/" << comm_world_rank() << ".csv"; 

  std::ofstream ofile;  
  ofile.open(stream.str().c_str());

  BOOST_FOREACH(Particle p, particles) {
    ofile << p[0] << ":" << p[1] << ",";
  } 

  ofile << 0 << ":" << 0;

  ofile.close();
#endif
/* an example from diy2-sfm
  diy::MemoryBuffer sendbuf;
  std::vector<int> sendcounts(np_group(), 0), sdispls(np_group(), 0);
  for (int i=0; i<np_group(); i++) {
    int s0 = sendbuf.size(); 
    diy::save(sendbuf, _finished_particles[i]);
    int s1 = sendbuf.size(); 
    sendcounts[i] = s1 - s0;
    sdispls[i] = s0;
  }

  diy::MemoryBuffer recvbuf; 
  std::vector<int> recvcounts(np_group(), 0), rdispls(np_group(), 0);

  MPI_Alltoall(
      (int*)sendcounts.data(), 
      1,
      MPI_INT, 
      (int*)recvcounts.data(), 
      1,
      MPI_INT, 
      _comm_group);

  int tot_size=0; 
  for (int i=0; i<np_group(); i++) {
    rdispls[i] = tot_size;
    tot_size += recvcounts[i];
  }

  recvbuf.buffer.resize(tot_size);

  MPI_Alltoallv(
      (void*)sendbuf.buffer.data(), 
      (int*)sendcounts.data(), 
      (int*)sdispls.data(),
      MPI_CHAR, 
      (void*)recvbuf.buffer.data(), 
      (int*)recvcounts.data(), 
      (int*)rdispls.data(), 
      MPI_CHAR, _comm_group);

  for (int i=0; i<np_group(); i++) {
    std::map<int, std::vector<Particle> > pts; // gid, pts
    diy::load(recvbuf, pts);

    // fprintf(stderr, "rank=%d, from rank %d, rdispl=%d, count=%d, mapsize=%d\n", 
    //     rank(), i, rdispls[i], recvcounts[i], pts.size());

    for (std::map<int, std::vector<Particle> >::iterator it = pts.begin(); it != pts.end(); it ++) {
      const int gid = it->first;
      Block &b = block(gid);
      for (std::vector<Particle>::iterator it1 = it->second.begin(); it1 != it->second.end(); it1 ++) {
        b.finished_particles[it1->id].push_back(*it1);
      }
    }
  }
*/
}

void CPTApp_Sync::add_workload()
{
  _round_steps ++;
  _epoch_steps ++;
  _total_steps ++;
}
