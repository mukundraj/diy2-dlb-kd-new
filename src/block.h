#ifndef _BLOCK_H
#define _BLOCK_H

#include <diy/master.hpp>
#include <diy/proxy.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/shared_mutex.hpp>
#include <boost/atomic.hpp>
#include "dlb.h"
#include "src/config.pb.h"

typedef diy::DiscreteBounds DBounds; 
typedef diy::ContinuousBounds CBounds;
typedef diy::RegularGridLink RGLink;

struct Block {
  Block();
  ~Block();

public:
  int gid;
  int lb[4], ub[4], // core lb/ub      start, end
      glb[4], gub[4]; // ghost lb/ub   start, end
  int lload[4], uload[4]; // bounds of data loading (added by Jiang)
  std::vector<float> nbr_bounds; // added by mraj
  std::vector<int> nbr_gids; // added by mraj
    std::vector<int> all_gids; // global gid list added by mraj
  std::vector<float> all_bounds; // global bounds list added by mraj
 
  // diy::RegularLink<CBounds> cons_kdtree_link; // added by mraj


  std::vector<float*> vars;
  //int neighbors[8]; // key=nbr_dir_t.  val is -1 if no neighbor

  /* input particles */
  //int num_particles;

  /* data bounds */
  DBounds data_bounds{4};
  CBounds core_bounds{4};

public:
  void get_ghost_st_sz(int ndims, int *gst, int *gsz);
  void get_load_st_sz(int ndims, int *lst, int *lsz);
  void get_ghost_load_st_sz(int ndims, int *gst, int *gsz, int *lst, int *lsz);
  //int neighbor_dir(int gid) const; // returns -1 if not neighbor
  void get_core_st_sz(int ndims, float *clb, float *cub);

public: // ``foreach'' funcs 
  void bil_add_block_2D(const PBDataset& d);
  void bil_add_block_3D(const PBDataset& d);
  void bil_add_block_4D(const PBDataset& d);

  void add_block_nek(const PBDataset& d);

  void show_geometry(const diy::Master::ProxyWithLink &cp, void*);
  void print_info();

public: // TOOD: ugly friend classes..
  friend class CPTApp;
  friend class CPTApp_Sync;
  friend class CPTApp_ThreadPool;
  friend class AddBlock;
  friend class VerifyBlockSync;
  friend class TraceBlockRoundRobin;
  friend class TraceBlockRoundRobin_Async;
  friend class TraceBlockAsync;
  friend class TraceBlockAsync_Interval;
  //friend class TraceBlockForKDTree;
  // stats
  size_t num_particles_initialized, num_particles_finished; 
  size_t num_particles_traced;
  bool initialized, finished;

  size_t num_particles_kd;
  //float *particles;
  std::vector<Particle> kd_particles;
  std::vector<Particle> kdtree_particles;
  std::vector<Particle> particles; // unfinished particles for k-d tree decomposition
  std::vector<Particle> pending_particles; // pending particles to be exchanged
  
  bool after_kdtree; 

  // workload
  boost::atomic<uint64_t> workload;  // number of steps.  total workload during the whole run

  // mutex
  boost::mutex mutex, mutex_init, mutex_finish;

};

#endif
