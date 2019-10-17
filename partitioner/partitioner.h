#ifndef _PARTITIONER_H
#define _PARTITIONER_H

// my own partitioner to probe maximum possible ghost size

#include <mpi.h>
#include <vector>
#include <stdint.h>


class RegularPartitioner
{
public: 
  RegularPartitioner(MPI_Comm comm, int ndims, int nruns, int nvars); 
  ~RegularPartitioner();

  void set_domain_size(const int*); 
  void set_num_blocks(int); 
  void set_given(const int*); 
  void set_ghost(int);
  void set_ghost(const int*);
  void set_ghost(const int*, const int*);

  void get_block_size(int*) const;
  void get_block_bounds(int *mins, int *maxs) const;

  void update_partition_pick_block(const int b, int* bst, int* bsz); // added by Jiang

  void update_partition();
  void update_ghosts(bool wrap);
  int64_t tot_block_mem() const;

  int global_num_blocks() const;
  int local_num_blocks() const;

  bool is_max_bounds() const;

public:
  int lid2gid(int lid) const; 
  int gid2lid(int gid) const; 
  int gid2rank(int gid) const; 
  int pt2gid(const float *pt) const; 
  // int pt2gid_core(const float *pt, std::vector<float>& nbr_bounds, std::vector<float>& nbr_gids) const; // added by mraj
  int pt2rank(const float *pt) const; 

public:
  void print_blk_info();

private: 
  const int _ndims, _nruns, _nvars;
  int _nblocks; 
  int _ds[4], 
      _given[4]; 
  int _lgs[4], _ugs[4];

  int _bs0[4]; // block size
  int _lat[4]; 

  MPI_Comm _comm; 
  int _rank, _np; 
 
private:
  typedef struct {
    int gid, lid; // global and local block ids
    int st[4], sz[4]; // starts and sizes
    int gst[4], gsz[4]; // starts and sizes with ghosts
  } block_t; 

  const block_t& block(int lid) const {return _blocks[lid];}
  block_t& block(int lid) {return _blocks[lid];}
  void block_apply_ghost(block_t& blk, const int *ds, const int *lgs, const int *ugs, bool wrap) const;

  std::vector<block_t> _blocks; // indexed by lid
};

#endif
