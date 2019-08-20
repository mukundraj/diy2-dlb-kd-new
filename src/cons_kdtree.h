#ifndef _CONS_KDTREE_H
#define _CONS_KDTREE_H

#include <diy/master.hpp>
#include <diy/proxy.hpp>
#include "dlb.h"
#include "block.h"

// constrained k-d tree decomposition with partial data duplication for SciVis'17 submission
// added by Jiang Zhang
// start from 12/06/2016

typedef diy::DiscreteBounds DBounds;
typedef diy::ContinuousBounds CBounds;
typedef diy::RegularContinuousLink RCLink; 
typedef diy::RegularGridLink RGLink;

struct ConstrainedKDTreeBlock
{
  struct Point // correspond to Particle in dlb.h
  {
    float&          operator[](unsigned i)              { return data[i]; }
    const float&    operator[](unsigned i) const        { return data[i]; }
    float           data[4];
    int             home_gid;
    int             num_steps;
    int             num_rounds;
    bool            finished;
  };
  std::vector<Point> points;
  //std::vector<Particle>  points;

  int domain_mins[4], domain_maxs[4]; 
};

struct PopulateMaster
{
  diy::Master* master;
  int          divs[4];
  int          dim;
  bool         space_only;
};

struct ExtractMaster
{
  diy::Master* master;
  int          dim;
  bool         space_only;
  int          divs[4];
  int          block_size[4];
  int          ghost_size[4];
  bool         constrained;
  bool         first;
  bool         async;
};

void populate_cons_kdtree_block(Block* blk, const diy::Master::ProxyWithLink& cp, void* aux_);
void extract_cons_kdtree_block(ConstrainedKDTreeBlock* b, const diy::Master::ProxyWithLink& cp, void* aux_);
double pt_cons_kdtree_exchange(
    diy::Master& master, const diy::Assigner& assigner, 
    std::vector<int> divisions, 
    const int dim, const bool space_only,
    int* block_size, int* ghost_size, 
    const bool constrained, const bool first, const bool async
  );

#endif
