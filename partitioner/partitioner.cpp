#include <iostream>
#include <cmath>
#include <cassert>
#include <boost/foreach.hpp>
#include "partitioner.h"

RegularPartitioner::RegularPartitioner(MPI_Comm comm, int ndims, int nruns, int nvars)
  : _comm(comm), 
    _ndims(ndims), _nruns(nruns), _nvars(nvars)
{
  memset(_given, 0, sizeof(int)*4); 
  MPI_Comm_size(comm, &_np); 
  MPI_Comm_rank(comm, &_rank); 
}

RegularPartitioner::~RegularPartitioner()
{
}

void RegularPartitioner::get_block_size(int *bs) const {
  memcpy(bs, _bs0, sizeof(int)*_ndims);
}

void RegularPartitioner::set_num_blocks(int nblocks)
{
  _nblocks = nblocks; 
}

void RegularPartitioner::set_domain_size(const int* ds)
{
  memcpy(_ds, ds, sizeof(int)*_ndims);
}

void RegularPartitioner::set_given(const int* given)
{
  memcpy(_given, given, sizeof(int)*_ndims);
}

void RegularPartitioner::set_ghost(int gs)
{
  for (int i=0; i<_ndims; i++) 
    _lgs[i] = _ugs[i] = gs;
}

void RegularPartitioner::set_ghost(const int *gs)
{
  set_ghost(gs, gs);
}

void RegularPartitioner::set_ghost(const int *lgs, const int *ugs)
{
  memcpy(_lgs, lgs, sizeof(int)*_ndims);
  memcpy(_ugs, ugs, sizeof(int)*_ndims);
}

void RegularPartitioner::update_partition_pick_block(const int b, int* bst, int* bsz) // added by Jiang
{
  int rem = _nblocks; 
  int bs0[4];  // expected block size
  int lat[4] = {0}; 
  int longest = 0; 

  for (int i=0; i<_ndims; i++) {
    if (_given[i] == 0) {
      lat[i] = 1; 
      bs0[i] = _ds[i];
    } else {
      lat[i] = _given[i]; 
      bs0[i] = _ds[i] / _given[i]; 
      assert(rem % _given[i] == 0); 
      rem /= _given[i]; 
    }
  } 

#if 1
  while (1) {
/*    // find the longest dimension
    longest = 0;
    for (int i=0; i<_ndims; i++) {
      if (_given[i] == 0 && bs0[i] > bs0[longest]) 
        longest = i; 
    }
*/
    // find smallest factor to subdivide the dimension
    for (int j=2; j<=rem; j++) {
      if (rem%j==0) {
        lat[longest] *= j; 
        bs0[longest] /= j; 
        rem /= j; 
        break;
      }
    } 
  
    longest = (longest+1)%_ndims;
    if (rem == 1) break; 
  }
#else // round-robin decomposition
  while (1) {
    // find smallest factor to subdivide the dimension
    for (int j=2; j<=rem; j++) {
      if (rem%j==0) {
        lat[longest] *= j; 
        bs0[longest] /= j; 
        rem /= j; 
        break;
      }
    } 
  
    if (rem == 1) break;

    // round-robin pick dimension
    longest = (longest+1)/_ndims;
  }
#endif
//fprintf(stderr, "lat=[%d, %d, %d, %d]\n", lat[0], lat[1], lat[2], lat[3]);
  // build blocks 
  int lid=0, gid=0; 
  int nl=_ndims<4?1:lat[3], // FIXME: only support to 4D
      nk=_ndims<3?1:lat[2], 
      nj=_ndims<2?1:lat[1], 
      ni=lat[0]; 
  int loc[4] = {0};
  int bs[4]; 

  loc[3]=0; 
  for (int l=0; l<nl; l++) 
  {
    loc[2]=0; 
    bs[3] = l==(nl-1) ? (_ds[3] - bs0[3]*(nl-1)) : bs0[3]; 
    for (int k=0; k<nk; k++) 
    {
      loc[1]=0; 
      bs[2] = k==(nk-1) ? (_ds[2] - bs0[2]*(nk-1)) : bs0[2]; 
      for (int j=0; j<nj; j++) 
      {
        loc[0]=0; 
        bs[1] = j==(nj-1) ? (_ds[1] - bs0[1]*(nj-1)) : bs0[1]; 
        for (int i=0; i<ni; i++) 
        {
          bs[0] = i==(ni-1) ? (_ds[0] - bs0[0]*(ni-1)) : bs0[0]; 

          if (gid == b) 
          { // pick one block
            block_t blk;
            blk.gid = gid; 
            blk.lid = lid ++; 
            for (int d=0; d<4; d++) { // MAX_DIM
              blk.st[d] = loc[d]; 
              blk.sz[d] = bs[d]; 

              bst[d] = loc[d];
              bsz[d] = bs[d];
            }
            _blocks.push_back(blk); 

            memcpy(_bs0, bs0, sizeof(int)*4); 
            memcpy(_lat, lat, sizeof(int)*4); 
            return;
          } 
          gid++; 
          loc[0] += bs0[0]; 
        }
        loc[1] += bs0[1]; 
      }
      loc[2] += bs0[2]; 
    }
    loc[3] += bs0[3]; 
  }
  
  memcpy(_bs0, bs0, sizeof(int)*4); 
  memcpy(_lat, lat, sizeof(int)*4);   
}

void RegularPartitioner::update_partition()
{
  int rem = _nblocks; 
  int bs0[4];  // expected block size
  int lat[4] = {0}; 
  int longest = 0; 

  for (int i=0; i<_ndims; i++) {
    if (_given[i] == 0) {
      lat[i] = 1; 
      bs0[i] = _ds[i];
    } else {
      lat[i] = _given[i]; 
      bs0[i] = _ds[i] / _given[i]; 
      assert(rem % _given[i] == 0); 
      rem /= _given[i]; 
    }
  } 

#if 1
  while (1) {
    // find the longest dimension
    longest = 0;
    for (int i=0; i<_ndims; i++) {
      if (_given[i] == 0 && bs0[i] > bs0[longest]) 
        longest = i; 
    }

    // find smallest factor to subdivide the dimension
    for (int j=2; j<=rem; j++) {
      if (rem%j==0) {
        lat[longest] *= j; 
        bs0[longest] /= j; 
        rem /= j; 
        break;
      }
    } 
  
    if (rem == 1) break; 
  }
#else // round-robin decomposition
  while (1) {
    // find smallest factor to subdivide the dimension
    for (int j=2; j<=rem; j++) {
      if (rem%j==0) {
        lat[longest] *= j; 
        bs0[longest] /= j; 
        rem /= j; 
        break;
      }
    } 
  
    if (rem == 1) break;

    // round-robin pick dimension
    longest = (longest+1)/_ndims;
  }
#endif

  // build blocks 
  int lid=0, gid=0; 
  int nl=_ndims<4?1:lat[3], // FIXME: only support to 4D
      nk=_ndims<3?1:lat[2], 
      nj=_ndims<2?1:lat[1], 
      ni=lat[0]; 
  int loc[4] = {0};
  int bs[4]; 

  loc[3]=0; 
  for (int l=0; l<nl; l++) 
  {
    loc[2]=0; 
    bs[3] = l==(nl-1) ? (_ds[3] - bs0[3]*(nl-1)) : bs0[3]; 
    for (int k=0; k<nk; k++) 
    {
      loc[1]=0; 
      bs[2] = k==(nk-1) ? (_ds[2] - bs0[2]*(nk-1)) : bs0[2]; 
      for (int j=0; j<nj; j++) 
      {
        loc[0]=0; 
        bs[1] = j==(nj-1) ? (_ds[1] - bs0[1]*(nj-1)) : bs0[1]; 
        for (int i=0; i<ni; i++) 
        {
          bs[0] = i==(ni-1) ? (_ds[0] - bs0[0]*(ni-1)) : bs0[0]; 

          if (gid%_np== _rank) 
          { // round robin
            block_t blk;
            blk.gid = gid; 
            blk.lid = lid ++; 
            for (int d=0; d<4; d++) { // MAX_DIM
              blk.st[d] = loc[d]; 
              blk.sz[d] = bs[d]; 
            }

            // block_print_info_par(blk, _rank); 
            _blocks.push_back(blk); 
          } 
          gid++; 
          loc[0] += bs0[0]; 
        }
        loc[1] += bs0[1]; 
      }
      loc[2] += bs0[2]; 
    }
    loc[3] += bs0[3]; 
  }
  
  memcpy(_bs0, bs0, sizeof(int)*4); 
  memcpy(_lat, lat, sizeof(int)*4); 
}

void RegularPartitioner::get_block_bounds(int *mins, int *maxs) const
{
  assert (_blocks.size() == 1);
  BOOST_FOREACH(const block_t& blk, _blocks) {
    for (int i = 0; i < _ndims; i ++) {
      mins[i] = blk.st[i];
      maxs[i] = blk.st[i] + blk.sz[i] - 1;
    }
  }
}

void RegularPartitioner::update_ghosts(bool wrap)
{
  BOOST_FOREACH(block_t& blk, _blocks)
    block_apply_ghost(blk, _ds, _lgs, _ugs, wrap); 
}

int RegularPartitioner::local_num_blocks() const
{
  return _blocks.size(); 
}

int RegularPartitioner::pt2rank(const float *pt) const
{
  int gid = pt2gid(pt); 
  return gid2rank(gid); 
}

int RegularPartitioner::lid2gid(int lid) const 
{
  return _blocks[lid].gid; 
}

int RegularPartitioner::gid2rank(int gid) const 
{
  return gid % _np; 
}

int RegularPartitioner::pt2gid(const float *pt) const 
{
  int gid; 
  int ids[_ndims]; 
  
  for (int i=0; i<_ndims; i++) 
    ids[i] = std::min((int)(pt[i]/_bs0[i]), (int)_lat[i]-1); // TODO: fail safe..

  switch (_ndims) 
  {
  case 2: gid = ids[1] * _lat[0] + ids[0]; break; 
  case 3: gid = ids[2] * _lat[1] * _lat[0] + ids[1] * _lat[0] + ids[0]; break; 
  case 4: gid = ids[3] * _lat[2] * _lat[1] * _lat[0] + ids[2] * _lat[1] * _lat[0] + ids[1] * _lat[0] + ids[0]; break;  
  default: assert(false); 
  }

  if (gid == -1) 
  {
    fprintf(stderr, "[pt2gid] error: pt={%f, %f, %f, %f}, bs={%d, %d, %d, %d}, indicies={%d, %d, %d, %d}\n", 
        pt[0], pt[1], pt[2], pt[3],
        _bs0[0], _bs0[1], _bs0[2], _bs0[3],
        ids[0], ids[1], ids[2], ids[3]); 
  }
  
  return gid; 
}

int RegularPartitioner::gid2lid(int gid) const 
{
  return gid / _np; 
}

bool RegularPartitioner::is_max_bounds() const
{
  BOOST_FOREACH(const block_t& blk, _blocks) {
    bool max = true;
    for (int i = 0; i < _ndims; i ++)
      max = max && blk.gst[i] == 0 && (blk.gst[i] + blk.gsz[i] == _ds[i]);
    if (max) return true;
  }

  return false;
}

void RegularPartitioner::block_apply_ghost(block_t& blk, const int *ds, const int *lgs, const int *ugs, bool wrap) const
{
#if 0
  int i; 
  for (i=0; i<_ndims; i++) {
    bool left  = blk.st[i] - lgs[i] >= 0,
         right = blk.st[i] + blk.sz[i] + ugs[i] <  ds[i];

    blk.gst[i] = left ? blk.st[i] - lgs[i] : blk.st[i];
    if (left && right) 
      blk.gsz[i] = blk.sz[i] + lgs[i] + ugs[i];
    else if (left && !right) 
      blk.gsz[i] = blk.sz[i] + lgs[i];
    else if (!left && right)
      blk.gsz[i] = blk.sz[i] + ugs[i];
    else 
      blk.gsz[i] = blk.sz[i];
  }
#else
  if (wrap) {
    for (int i = 0; i < _ndims; i ++) {
      blk.gst[i] = blk.st[i] - lgs[i];
      blk.gsz[i] = blk.st[i] + blk.sz[i] - blk.gst[i] + ugs[i];
    }
  } else {
    for (int i = 0; i < _ndims; i ++) {
      blk.gst[i] = std::max(0, blk.st[i] - lgs[i]);
      blk.gsz[i] = std::min(ds[i] - blk.gst[i], blk.st[i] + blk.sz[i] - blk.gst[i] + ugs[i]);
    }
  }
#endif
}

int64_t RegularPartitioner::tot_block_mem() const
{
  int64_t mem = 0;
  BOOST_FOREACH(const block_t& blk, _blocks) {
    int64_t count = 1;
    for (int i=0; i<_ndims; i++) {
      count *= blk.gsz[i];
    }
    count *= _nvars * _nruns;
    mem += count * sizeof(float);
    // fprintf(stderr, "gsz={%d, %d, %d, %d}, nv=%d, %lld\n", blk.gsz[0], blk.gsz[1], blk.gsz[2], blk.gsz[3], _nvars, mem);
  }
  return mem;
}

void RegularPartitioner::print_blk_info()
{
  BOOST_FOREACH(const block_t& blk, _blocks) {
    fprintf(stderr, "rank=%d, start={%d, %d, %d, %d}, end={%d, %d, %d, %d}, size={%d, %d, %d, %d}, gstart={%d, %d, %d, %d}, gend={%d, %d, %d, %d}\n", 
        _rank, blk.st[0], blk.st[1], blk.st[2], blk.st[3], blk.st[0]+blk.sz[0]-1, blk.st[1]+blk.sz[1]-1, blk.st[2]+blk.sz[2]-1, blk.st[3]+blk.sz[3]-1,
        blk.sz[0], blk.sz[1], blk.sz[2], blk.sz[3], 
        blk.gst[0], blk.gst[1], blk.gst[2], blk.gst[3], blk.gst[0]+blk.gsz[0]-1, blk.gst[1]+blk.gsz[1]-1, blk.gst[2]+blk.gsz[2]-1, blk.gst[3]+blk.gsz[3]-1);
  }
}
