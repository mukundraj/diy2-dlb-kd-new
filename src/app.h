#ifndef _CPTAPP_H
#define _CPTAPP_H

#include <diy/mpi.hpp>
#include <diy/master.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/decomposition.hpp>
#include <diy/mpi/datatypes.hpp>
#include "dlb.h"
#include "block.h"
#include "src/config.pb.h"

class CPTApp {
  typedef diy::DiscreteBounds DBounds;
  typedef diy::ContinuousBounds CBounds;
  typedef diy::RegularGridLink RGLink;

public:
  CPTApp();
  ~CPTApp();

  virtual void init(int argc, char **argv);
  virtual void deinit();
  virtual void exec() = 0;

public:
  MPI_Comm _comm_world; 
  int _comm_world_rank, _comm_world_size;

  int comm_world_rank() const { return _comm_world_rank; }
  int comm_world_size() const { return _comm_world_size; }
  int rank(int gid) const { return _assigner->rank(gid); }

  MPI_Comm comm_world() const { return _comm_world; }

  void init_mpi();

public:
  int num_dims() const { return _num_dims; }
  int num_threads() const { return _num_threads; }
  int num_blocks() const { return _num_blocks; }

  const int* domain_size() const { return _domain_size; }

  Block& block(int gid) {
    int lid = _master->lid(gid); 
    return *static_cast<Block*>(_master->block(lid));
  }
  const std::vector<int>& gids() const {return _gids;}

  int lid2gid(int lid) const {return _gids[lid];}
  int gid2lid(int gid) const {return _master->lid(gid);}

  int idx2id(const int idx[4]) const;

  int pt2gid(const float *X) const {return _decomposer->point_to_gid(X);} // TODO 
  bool inside_domain(const float *X) const {
    for (int i=0; i<_num_dims; i++) 
      if (X[i]<0 || X[i]>=_domain_size[i]-1) return false;
    return true;
  }

  int bound_gid(const int id) const { 
    if (id >= _bound_gids.size()) assert(false); 
    return _bound_gids[id]; 
  }
  bool is_kd_tree() const { return _appconf.kd_tree(); }
  bool space_only() const { 
    if (num_dims() <= 3) return false;
    else return _appconf.space_only();
  }

  int pred_val() const { return _appconf.prediction(); }

  bool is_ptinblock(const Block& b, const float* pt);

protected:
  bool parse_arguments(int argc, char **argv);
  bool parse_config(const std::string& filename);
  void bcast_config();

  void probe_bs(int *bsz) const;
  void probe_elastic_ghosts(int *bsz, int *ghosts) const; 
  void probe_elastic_ghost_block_size(const int factor, int *bsz, int *gsz) const;
  void regular_probe_elastic_ghost_size(int *gs) const;

  PBConfig _appconf;
  PBDataset _dataconf;

protected: // for final statistics
  int max_num(const std::vector<int>& nums);
  int min_num(const std::vector<int>& nums);

  float avg_num(const std::vector<int>& nums);
  float var_num(const std::vector<int>& nums);

  int sum_num(const std::vector<int>& nums);

protected:
  virtual void stat();
  void write_output_file();

  double _timestamp_start, _timestamp_init, _timestamp_finish, _time_io;

protected:
  std::vector<double> _timestamps;
  std::vector<int> _timecategories;

  std::vector<float> _balance;

  bool _constrained;

protected:
  int _num_dims, _num_threads;
  int _num_blocks, _nb_per_proc; 
  int _domain_size[4], _ghost_size[4];
  int _block_size[4];

  diy::mpi::communicator *_communicator;
  diy::Master *_master;
  diy::RoundRobinAssigner *_assigner;
  diy::RegularDecomposer<DBounds> *_decomposer;
  DBounds _domain;

  int _nt_loops;
  std::vector<int> _divisions;
  std::vector<int> _bound_gids;

  std::vector<int> _gids; // lid2gid

};

#endif
