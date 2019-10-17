#ifndef DIY_DETAIL_ALGORITHMS_CONSTRAINED_KDTREE_HPP
#define DIY_DETAIL_ALGORITHMS_CONSTRAINED_KDTREE_HPP

#include <vector>
#include <cassert>
#include <algorithm>
#include "../../partners/all-reduce.hpp"

// constrained k-d tree decomposition with partial data duplication for SciVis'17 submission
// added by Jiang Zhang
// start from 12/06/2016

#define DEBUG 0

namespace diy
{
namespace detail
{

struct ConstrainedKDTreePartners;

template<class Block, class Point>
struct ConstrainedKDTreePartition
{
    typedef     diy::RegularContinuousLink      RCLink;
    typedef     diy::ContinuousBounds           Bounds;

    typedef     std::vector<size_t>             Histogram;

                ConstrainedKDTreePartition(int                             dim,
                                std::vector<Point>  Block::*    points,
                                const int*                      block_size, // added by Jiang
                                const int*                      ghost_size, // added by Jiang                               
                                const bool                      constrained, // added by Jiang
                                const bool                      first,      // added by Jiang
                                //const std::vector<int>          divisions,  // added by Jiang
                                size_t                          bins):
                    dim_(dim), 
                    constrained_(constrained),
                    first_(first),
                    //divisions_(divisions),
                    points_(points), 
                    bins_(bins)            {
                        for (int i = 0; i < dim; i ++) { // added by Jiang
                            block_size_[i] = block_size[i];
                            ghost_size_[i] = ghost_size[i];
                        }

                        prev_dim_ = 0;
                    }

    void        operator()(void* b_, const diy::ReduceProxy& srp, const ConstrainedKDTreePartners& partners) const;

    int         divide_gid(int gid, bool lower, int round, int rounds) const;
    void        update_links(Block* b, const diy::ReduceProxy& srp, int dim, int round, int rounds, const Bounds& domain) const;
    void        split_to_neighbors(Block* b, const diy::ReduceProxy& srp, int dim) const;
    //diy::Direction
    //            find_wrap(const Bounds& bounds, const Bounds& nbr_bounds, const Bounds& domain) const;

    void        compute_local_histogram(Block* b, const diy::ReduceProxy& srp, int dim) const;
    void        add_histogram(Block* b, const diy::ReduceProxy& srp, Histogram& histogram) const;
    void        receive_histogram(Block* b, const diy::ReduceProxy& srp,       Histogram& histogram) const;
    void        forward_histogram(Block* b, const diy::ReduceProxy& srp, const Histogram& histogram) const;

    void        enqueue_exchange(Block* b, const diy::ReduceProxy& srp, int dim, const Histogram& histogram) const;
    void        dequeue_exchange(Block* b, const diy::ReduceProxy& srp, int dim) const;

    void        update_neighbor_bounds(Bounds& bounds, float split, int dim, bool lower) const;
    bool        intersects(const Bounds& x, const Bounds& y, int dim, const Bounds& domain) const;
    float       find_split(const Bounds& changed, const Bounds& original) const;

    //float       find_middle(const int& dim, const int& block_size, const Bounds& bbox) const; // compute middle of bbox

/*    
    float       find_middle(
                    const int& block_size, 
                    const int* bbox_minmax
                  ) const;

    void        fill_bounds(
                    const int& ghost_size, 
                    const int* domain_minmax, 
                    const int* core, 
                    int* bounds
                  ) const;
*/

    void        find_middle_and_range(
                    const int rank,
                    const int& dim,
                    const float* divs,
                    const int* domain_minmax, 
                    //const float* bbox_minmax, 
                    const float* cutbox_minmax,
                    const int& block_size,
                    const int& ghost_size,
                    float& mid_div,
                    //float& middle, 
                    float* range
                  ) const;

    int                             dim_;
    std::vector<Point>  Block::*    points_;
    int                             block_size_[4];
    int                             ghost_size_[4];
    size_t                          bins_;
    bool                            first_;
    bool                            constrained_;
    //std::vector<int>                divisions_;
    int                             prev_dim_;  // added by Jiang
};

}
}

struct diy::detail::ConstrainedKDTreePartners
{
  // bool = are we in a swap (vs histogram) round
  // int  = round within that partner
  typedef           std::pair<bool, int>                    RoundType;
  typedef           diy::ContinuousBounds                   Bounds;

                    ConstrainedKDTreePartners(int dim, int nblocks, const Bounds& domain_, const std::vector<int> divisions_):
                        decomposer(1, interval(0,nblocks-1), nblocks),
                        histogram(decomposer, 2),
                        swap(decomposer, 2, false),
                        domain(domain_),
                        divisions(divisions_)
  {
#if 1
    int rdim = 0;
    std::vector<int> divs = divisions_;

    for (unsigned i = 0; i < swap.rounds(); ++i)
    {
      while (divs[rdim] / 2 < 1) rdim = (rdim + 1) % dim; 
      divs[rdim] /= 2;

      // fill histogram rounds (histogram for reduce)
      for (unsigned j = 0; j < histogram.rounds(); ++j)
      {
        rounds_.push_back(std::make_pair(false, j));
        //dim_.push_back(i % dim);
        dim_.push_back(rdim % dim); // i % dim
        if (j == histogram.rounds() / 2 - 1 - i)
            j += 2*i;
      }

      // fill swap round 
      rounds_.push_back(std::make_pair(true, i));
      //dim_.push_back(i % dim);
      dim_.push_back(rdim % dim); // i % dim

      // fill link round
      rounds_.push_back(std::make_pair(true, -1));          // (true, -1) signals link round
      //dim_.push_back(i % dim);
      dim_.push_back(rdim % dim); // i % dim

      rdim = (rdim + 1) % dim;
    }
#else
    for (unsigned i = 0; i < swap.rounds(); ++i)
    {
      // fill histogram rounds
      for (unsigned j = 0; j < histogram.rounds(); ++j)
      {
        rounds_.push_back(std::make_pair(false, j));
        dim_.push_back(i % dim);
        if (j == histogram.rounds() / 2 - 1 - i)
            j += 2*i;
      }

      // fill swap round
      rounds_.push_back(std::make_pair(true, i));
      dim_.push_back(i % dim);

      // fill link round
      rounds_.push_back(std::make_pair(true, -1));          // (true, -1) signals link round
      dim_.push_back(i % dim);
    }
#endif
  }

  size_t           rounds() const                              { return rounds_.size(); }
  size_t           swap_rounds() const                         { return swap.rounds(); }

  std::vector<int> get_divisions() const                       { return divisions; }

  int              dim(int round) const                        { return dim_[round]; }
  bool             swap_round(int round) const                 { return rounds_[round].first; }
  int              sub_round(int round) const                  { return rounds_[round].second; }

  inline bool      active(int round, int gid, const diy::Master& m) const
  {
    if (round == (int) rounds())      // final round ?
        return true;
    else if (swap_round(round) && sub_round(round) < 0)     // link round
        return true;
    else if (swap_round(round))     // swap round
        return swap.active(sub_round(round), gid, m);
    else      // histogram round
        return histogram.active(sub_round(round), gid, m);
  }

  inline void      incoming(int round, int gid, std::vector<int>& partners, const diy::Master& m) const
  {
    if (round == (int) rounds()) {       // final round ?
        link_neighbors(-1, gid, partners, m);
    }
    else if (swap_round(round) && sub_round(round) < 0) {       // link round (before dequeue_exchange)
        swap.incoming(sub_round(round - 1) + 1, gid, partners, m);
    }
    else if (swap_round(round)) {      // swap round (before receive_histogram and enqueue_exchange)
        histogram.incoming(histogram.rounds(), gid, partners, m);
    }
    else        // histogram round
    {
        if (round > 0 && sub_round(round) == 0) {    // first histogram round during a new swap round (before compute local histogram)
            link_neighbors(-1, gid, partners, m);
        }
        else if (round > 0 && sub_round(round - 1) != sub_round(round) - 1) {       // jump through the histogram rounds (other times)
            histogram.incoming(sub_round(round - 1) + 1, gid, partners, m);
        }
        else {   // before add histogram (the first time)
            histogram.incoming(sub_round(round), gid, partners, m);
        }
    }
  }

  inline void      outgoing(int round, int gid, std::vector<int>& partners, const diy::Master& m) const
  {
    if (round == (int) rounds()) {      // final round
        swap.outgoing(sub_round(round-1) + 1, gid, partners, m);
    }
    else if (swap_round(round) && sub_round(round) < 0) {      // link round (before dequeue_exchange)
        link_neighbors(-1, gid, partners, m);
    }
    else if (swap_round(round)) {       // swap round (before receive_histogram and enqueue_exchange)
        swap.outgoing(sub_round(round), gid, partners, m);
    }
    else {       // histogram round (before compute local histogram, add histogram)
        histogram.outgoing(sub_round(round), gid, partners, m);
    }
  }

  inline void      link_neighbors(int, int gid, std::vector<int>& partners, const diy::Master& m) const
  {
    int         lid  = m.lid(gid);
    diy::Link*  link = m.link(lid);

    std::set<int> result;       // partners must be unique
    for (int i = 0; i < link->size(); ++i)
        result.insert(link->target(i).gid);

    for (std::set<int>::const_iterator it = result.begin(); it != result.end(); ++it)
        partners.push_back(*it);
  }

  // 1-D domain to feed into histogram and swap
  diy::RegularDecomposer<diy::DiscreteBounds>   decomposer;

  diy::RegularAllReducePartners     histogram;
  diy::RegularSwapPartners          swap;

  std::vector<RoundType>            rounds_;
  std::vector<int>                  dim_;

  Bounds                            domain;
  std::vector<int>                  divisions;
};

template<class Block, class Point>
void
diy::detail::ConstrainedKDTreePartition<Block,Point>::
operator()(void* b_, const diy::ReduceProxy& srp, const ConstrainedKDTreePartners& partners) const
{
    Block* b = static_cast<Block*>(b_);



    int dim;
    if (srp.round() < partners.rounds()) {
// before the final
        dim = partners.dim(srp.round());
    }
    else {
// only in the final
        dim = partners.dim(srp.round() - 1);
    }

    if (srp.round() == partners.rounds()) {
// only in the final 
        // this link update is used for local communication after decomposition TODO
        update_links(b, srp, dim, partners.sub_round(srp.round() - 2), partners.swap_rounds(), partners.domain); // -1 would be the "uninformative" link round
    }
    else if (partners.swap_round(srp.round()) && partners.sub_round(srp.round()) < 0)       // link round
    {
// after each time operator 3 (enqueue_exchange)
        dequeue_exchange(b, srp, dim);         // from the swap round
        split_to_neighbors(b, srp, dim);    
    }
    else if (partners.swap_round(srp.round()))        // swap round
    {
// after each time operator 6 (some procs forward_histogram)
        Histogram   histogram;
        receive_histogram(b, srp, histogram);
        enqueue_exchange(b, srp, dim, histogram);    // swap
    } else if (partners.sub_round(srp.round()) == 0)       // histogram rounds below
    {
        if (srp.round() > 0) // after each iteration of exchanging particles (dequeue exchange)
        {
/*
            int prev_dim = dim - 1;
            if (prev_dim < 0)
                prev_dim += dim_;
            update_links(b, srp, prev_dim, partners.sub_round(srp.round() - 2), partners.swap_rounds(), partners.domain);    // -1 would be the "uninformative" link round
*/
            int prev_dim = partners.dim(srp.round() - 1);
            // this link update is used for local communication after decomposition TODO
            update_links(b, srp, prev_dim, partners.sub_round(srp.round() - 2), partners.swap_rounds(), partners.domain);    // -1 would be the "uninformative" link round            
        }
// compute local histogram each time in new cut box 
        compute_local_histogram(b, srp, dim);
    } else if (partners.sub_round(srp.round()) < (int) partners.histogram.rounds()/2)    // lower half of the group 
    {
        Histogram   histogram(bins_);
        add_histogram(b, srp, histogram);
        srp.enqueue(srp.out_link().target(0), histogram);
    }
    else     // higher half of the group
    {
        Histogram   histogram(bins_);
        add_histogram(b, srp, histogram);
        forward_histogram(b, srp, histogram);
    }

    // fprintf(stderr, "DB round %d gid %d size %ld\n", srp.round(), srp.gid(), b->points.size());
}

template<class Block, class Point>
int
diy::detail::ConstrainedKDTreePartition<Block,Point>::
divide_gid(int gid, bool lower, int round, int rounds) const
{
    if (lower)
        gid &= ~(1 << (rounds - 1 - round));
    else
        gid |=  (1 << (rounds - 1 - round));
    return gid;
}

// round here is the outer iteration of the algorithm
template<class Block, class Point>
void
diy::detail::ConstrainedKDTreePartition<Block,Point>::
update_links(Block* b, const diy::ReduceProxy& srp, int dim, int round, int rounds, const Bounds& domain) const
{
    int         gid  = srp.gid();
    int         lid  = srp.master()->lid(gid);
    RCLink*     link = static_cast<RCLink*>(srp.master()->link(lid));

    // (gid, dir) -> i
    std::map<std::pair<int,diy::Direction>, int> link_map;
    for (int i = 0; i < link->size(); ++i)
        link_map[std::make_pair(link->target(i).gid, link->direction(i))] = i;

    // NB: srp.enqueue(..., ...) should match the link
    std::vector<float>  splits(link->size());
    for (int i = 0; i < link->size(); ++i)
    {
        float split; diy::Direction dir(dim_, 0);

        int in_gid = link->target(i).gid;
        while(srp.incoming(in_gid))
        {
            srp.dequeue(in_gid, split);
            srp.dequeue(in_gid, dir);

            // reverse dir
            for (int j = 0; j < dim_; ++j)
                dir[j] = -dir[j];

            int k = link_map[std::make_pair(in_gid, dir)];
            splits[k] = split;
        }
    }

    RCLink      new_link(dim_, link->core(), link->core());
    //new_link.set_regular_core(link->regular_core());
    new_link.set_divs(link->divs());

    bool lower = !(gid & (1 << (rounds - 1 - round)));

    // fill out the new link
    for (int i = 0; i < link->size(); ++i)
    {
        diy::Direction  dir      = link->direction(i);
        //diy::Direction  wrap_dir = link->wrap(i);     // we don't use existing wrap, but restore it from scratch
        if (dir[dim] != 0)
        {
            if ((dir[dim] < 0 && lower) || (dir[dim] > 0 && !lower))
            {
                int nbr_gid = divide_gid(link->target(i).gid, !lower, round, rounds);
                diy::BlockID nbr = { nbr_gid, srp.assigner().rank(nbr_gid) };
                new_link.add_neighbor(nbr);

                new_link.add_direction(dir);

                Bounds bounds = link->bounds(i);
                update_neighbor_bounds(bounds, splits[i], dim, !lower);
                new_link.add_bounds(bounds);

                //if (wrap)
                //    new_link.add_wrap(find_wrap(new_link.bounds(), bounds, domain));
                //else
                new_link.add_wrap(diy::Direction(dim_, 0));
            }
        } else // non-aligned side
        {
            for (int j = 0; j < 2; ++j)
            {
                int nbr_gid = divide_gid(link->target(i).gid, j == 0, round, rounds);

                Bounds  bounds  = link->bounds(i);
                update_neighbor_bounds(bounds, splits[i], dim, j == 0);

                if (intersects(bounds, new_link.bounds(), dim, domain))
                {
                    diy::BlockID nbr = { nbr_gid, srp.assigner().rank(nbr_gid) };
                    new_link.add_neighbor(nbr);
                    new_link.add_direction(dir);
                    new_link.add_bounds(bounds);

                    //if (wrap)
                    //    new_link.add_wrap(find_wrap(new_link.bounds(), bounds, domain));
                    //else
                    new_link.add_wrap(diy::Direction(dim_, 0));
                }
            }
        }
    }

    // add link to the dual block
    int dual_gid = divide_gid(gid, !lower, round, rounds);
    diy::BlockID dual = { dual_gid, srp.assigner().rank(dual_gid) };
    new_link.add_neighbor(dual);

    Bounds nbr_bounds = link->bounds();     // old block bounds
    update_neighbor_bounds(nbr_bounds, find_split(new_link.bounds(), nbr_bounds), dim, !lower);
    new_link.add_bounds(nbr_bounds);

    new_link.add_wrap(diy::Direction());    // dual block cannot be wrapped

    if (lower)
    {
        diy::Direction right(dim_, 0);
        right[dim] = 1;
        new_link.add_direction(right);
    } else
    {
        diy::Direction left(dim_, 0);
        left[dim] = -1;
        new_link.add_direction(left);
    }
    // update the link; notice that this won't conflict with anything since
    // reduce is using its own notion of the link constructed through the
    // partners
    link->swap(new_link);   
}

template<class Block, class Point>
void
diy::detail::ConstrainedKDTreePartition<Block,Point>::
split_to_neighbors(Block* b, const diy::ReduceProxy& srp, int dim) const
{
    int         lid  = srp.master()->lid(srp.gid());
    RCLink*     link = static_cast<RCLink*>(srp.master()->link(lid));
    float split = find_split(link->core(), link->bounds());

    for (int i = 0; i < link->size(); ++i)
    {
        srp.enqueue(link->target(i), split);
        srp.enqueue(link->target(i), link->direction(i));
    }
}

template<class Block, class Point>
void
diy::detail::ConstrainedKDTreePartition<Block,Point>::
compute_local_histogram(Block* b, const diy::ReduceProxy& srp, int dim) const
{
    int         lid  = srp.master()->lid(srp.gid());
    RCLink*     link = static_cast<RCLink*>(srp.master()->link(lid));

    // compute and enqueue local histogram
    Histogram histogram(bins_);

    // link core is the cutting box (initially is the entire domain)
    float   width = (link->core().max[dim] - link->core().min[dim])/bins_; 
    for (size_t i = 0; i < (b->*points_).size(); ++i)
    {
        float x = (b->*points_)[i][dim];
        int loc = (x - link->core().min[dim]) / width;
        if (loc < 0)
        {
            std::cerr << loc << " " << x << " " << link->core().min[dim] << std::endl;
            std::abort();
        }
        if (loc >= (int) bins_)
            loc = bins_ - 1;
        // ++(histogram[loc]);
        int wgt = (b->*points_)[i].wgt; // mraj
        (histogram[loc])+=wgt;
        // fprintf(stderr, "%d ", wgt);
    }

    srp.enqueue(srp.out_link().target(0), histogram);
}

template<class Block, class Point>
void
diy::detail::ConstrainedKDTreePartition<Block,Point>::
add_histogram(Block* b, const diy::ReduceProxy& srp, Histogram& histogram) const
{
    // dequeue and add up the histograms
    for (int i = 0; i < srp.in_link().size(); ++i)
    {
        int nbr_gid = srp.in_link().target(i).gid;

        Histogram hist;
        srp.dequeue(nbr_gid, hist);
        for (size_t i = 0; i < hist.size(); ++i)
            histogram[i] += hist[i];
    }
}

template<class Block, class Point>
void
diy::detail::ConstrainedKDTreePartition<Block,Point>::
receive_histogram(Block* b, const diy::ReduceProxy& srp, Histogram& histogram) const
{
    srp.dequeue(srp.in_link().target(0).gid, histogram);
}

template<class Block, class Point>
void
diy::detail::ConstrainedKDTreePartition<Block,Point>::
forward_histogram(Block* b, const diy::ReduceProxy& srp, const Histogram& histogram) const
{
    for (int i = 0; i < srp.out_link().size(); ++i) {
        srp.enqueue(srp.out_link().target(i), histogram);
    }
}
/*
template<class Block, class Point>
float
diy::detail::ConstrainedKDTreePartition<Block,Point>::
find_middle(
    const int& block_size, 
    const int* bbox_minmax
  ) const
{
#if 0
  int middle = std::ceil((float)(bbox.max[dim] + bbox.min[dim] ) / 2.);

  return (float)middle;
#else
  float bound = bbox_minmax[1] - bbox_minmax[0] + 1;

  float middle = bbox_minmax[0];
  int size = std::floor(bound / (float)block_size);
  if (!(size >= 2 && size % 2 == 0)) {
    fprintf(stderr, "bbox bound= {%d, %d}, block size = %d, size = %d\n", 
        bbox_minmax[0], bbox_minmax[1], block_size, size);
    assert (false);
  }
  for (int i = 0; i < size/2; i ++)
    middle += block_size;
  
  return middle;
#endif
}
*/

/*
template<class Block, class Point>
void
diy::detail::ConstrainedKDTreePartition<Block,Point>::
fill_bounds(
    const int& ghost_block_size, 
    const int* domain_minmax, 
    const int* core, 
    int* bounds
  ) const
{
  float bounds_max, bounds_min;
  float add = (float)ghost_block_size/2. - (float)(core[1]-core[0]+1)/2.;
  assert (add >= 0);
  bounds_max = core[1] + add;
  bounds_min = core[0] - add;

  if (ghost_block_size > (domain_minmax[1] - domain_minmax[0] + 1)) {
    bounds_max = domain_minmax[1];
    bounds_min = domain_minmax[0];
  } else {
    if (bounds_max > domain_minmax[1]) {
      bounds_max = domain_minmax[1];
      bounds_min = bounds_max - ghost_block_size + 1;
    }

    if (bounds_min < domain_minmax[0]) {
      bounds_min = domain_minmax[0];
      bounds_max = bounds_min + ghost_block_size - 1;
    }
  }

  bounds[1] = std::floor(bounds_max);
  bounds[0] = std::floor(bounds_min);
  assert (bounds[1] >= core[1] && bounds[0] <= core[0]);
}
*/
template<class Block, class Point>
void 
diy::detail::ConstrainedKDTreePartition<Block,Point>::
find_middle_and_range(
    const int rank,
    const int& dim,
    const float* divs,
    const int* domain_minmax, 
    //const float* bbox_minmax, 
    const float* cutbox_minmax,
    const int& block_size,
    const int& ghost_size,
    float& mid_div,
    //float& middle, 
    float* range
  ) const
{
#if 0
  float bound = bbox_minmax[1] - bbox_minmax[0] + 1;

  middle = bbox_minmax[0];
  int size = std::floor(bound / (float)block_size);
  if (!(size >= 2 && size % 2 == 0)) {
    if (bbox_minmax[1] == domain_minmax[1] && (size-1) % 2 == 0) {
      size = size - 1;
    } else {
      fprintf(stderr, "dim = %d, bbox bound = {%f, %f}, block size = %d, size = %d\n", 
        dim, bbox_minmax[0], bbox_minmax[1], block_size, size);
      assert (false);
    }
  }
#else
  int size = std::floor(divs[1] - divs[0] + 1);
  if (!(size >= 2 && size % 2 == 0)) 
  {
    fprintf(stderr, "dim = %d, divs = {%f, %f}, block size = %d, size = %d\n", 
        dim, divs[0], divs[1], block_size, size);
    assert (false);
  }
#endif
  //for (int i = 0; i < size/2; i ++) middle += block_size;
  //middle = bbox_minmax[0];
  int bbox_min = divs[0] * block_size;
  mid_div = divs[0];

  int distance = size/2;
  //middle += block_size * distance;
  mid_div += (float)distance;
   
  int left[2], right[2], ghost_left[2], ghost_right[2];
  int left_min = (distance + std::floor(divs[0]) - 1) * block_size;
  left[0] = std::max(left_min, bbox_min); //(int)bbox_minmax[0]);
  left[1] = left[0] + block_size - 1;
  right[0] = left[1] + 1;
  right[1] = std::min(right[0] + block_size - 1, domain_minmax[1]);

  //fill_bounds(ghost_block_size, domain_minmax, left, ghost_left);
  //fill_bounds(ghost_block_size, domain_minmax, right, ghost_right);
  //ghost_left[0] = std::max(domain_minmax[0], left[0] - ghost_size[i]);
  ghost_left[1] = std::min(domain_minmax[1], left[1] + ghost_size);

  ghost_right[0] = std::max(domain_minmax[0], right[0] - ghost_size);
  //ghost_right[1] = std::min(domain_minmax[1], right[1] + ghost_size[i]);

  range[0] = std::max(cutbox_minmax[0], (float)ghost_right[0]);
  range[1] = std::min(cutbox_minmax[1], (float)ghost_left[1]);

  if (range[0] > range[1] || range[0] < cutbox_minmax[0] || range[1] > cutbox_minmax[1]) {
    fprintf(stderr, "Errror: range [%f, %f] exceeds cutbox [%f, %f]\n", range[0], range[1],
        cutbox_minmax[0], cutbox_minmax[1]);
    assert (false);
  }
}

// add constrains according to the ghost layers (added by Jiang Zhang)
template<class Block, class Point>
void
diy::detail::ConstrainedKDTreePartition<Block,Point>::
enqueue_exchange(Block* b, const diy::ReduceProxy& srp, int dim, const Histogram& histogram) const 
{//fprintf(stderr, "start to enqueue_exchange.\n");
    int         lid  = srp.master()->lid(srp.gid());
    RCLink*     link = static_cast<RCLink*>(srp.master()->link(lid));

    int rank = srp.master()->communicator().rank();

    int k = srp.out_link().size();

    if (k == 0)        // final round; nothing needs to be sent; this is actually redundant
        return;

    // find the split to split points to two parts with equal number 
    //float middle = 0; 
    float split = 0, mid_div = 0;
    if (constrained_) { // constrained k-d tree decomposition
      // link->regular_core(): regular bounding box (bbox)
      // link->core(): split bounding box (cut_bbox)
      // determine range of split

      int domain_minmax[2] = {b->domain_mins[dim], b->domain_maxs[dim]}; 
      float divs[2] = {link->divs().min[dim], link->divs().max[dim]};
      //float bbox_minmax[2] = {link->regular_core().min[dim], link->regular_core().max[dim]},
      float cutbox_minmax[2] = {link->core().min[dim], link->core().max[dim]};
      float range[2];
      //find_middle_and_range(dim, b->domain_mins[dim], b->domain_maxs[dim], link->regular_core(), middle, range);

      find_middle_and_range(
          rank,
          dim,
          divs,
          domain_minmax, 
          //bbox_minmax, 
          cutbox_minmax,
          block_size_[dim],
          ghost_size_[dim],
          mid_div,
          //middle, 
          range
        ); 

      if (first_) split = mid_div * block_size_[dim];
      else {
        size_t total = 0;
        for (size_t i = 0; i < histogram.size(); ++i)
          total += histogram[i];

        size_t cur   = 0;
        float  width = (link->core().max[dim] - link->core().min[dim])/bins_;
        size_t front = 0;
        while (link->core().min[dim] + width*front < range[0]) {
          cur += histogram[front]; 
          front ++;
        }

        //float  split = 0;
        split = range[0];
        if (cur < total/2) 
        {
          for (size_t i = front; i < histogram.size(); i ++) {
            if (cur + histogram[i] > total/2) {
              split = link->core().min[dim] + width*i;
              break;
            }
            cur += histogram[i];
          }
          split = std::min(range[1], split); 
        }

        if (split == link->core().min[dim] || split == link->core().max[dim]) {
          //split = (link->core().min[dim] + link->core().max[dim]) / 2.f;
          //fprintf(stderr, "Error: dim %d, split %f exceeds range [%f, %f], core [%f, %f]\n", 
          //    dim, split, range[0], range[1], link->core().min[dim], link->core().max[dim]);
          //assert (false);
          split = (range[0] + range[1]) / 2.f;
        }
/*
        if (split > range[1] || split < range[0]) {
          fprintf(stderr, "Error: split %f exceeds range [%f, %f], core [%f, %f]\n", 
              split, range[0], range[1], link->core().min[dim], link->core().max[dim]);
          assert (false);
        }
*/
        //fprintf(stderr, "rank %d, dim=%d, split=%f, range=[%f,%f], core=[%f,%f]\n", 
        //  rank, dim, split, range[0], range[1], link->core().min[dim], link->core().max[dim]);
      }
    } else { // regular k-d tree decomposition
      // pick split points
      size_t total = 0;
      for (size_t i = 0; i < histogram.size(); ++i)
        total += histogram[i];
      //fprintf(stderr, "Histogram total: %lu\n", total);

      size_t cur   = 0;
      float  width = (link->core().max[dim] - link->core().min[dim])/bins_;
      //float  split = 0;
      for (size_t i = 0; i < histogram.size(); ++i) {
        if (cur + histogram[i] > total/2) {
            split = link->core().min[dim] + width*i;
            break;
        }
        cur += histogram[i];
      }

      if (split == link->core().min[dim] || split == link->core().max[dim])
        split = (link->core().min[dim] + link->core().max[dim]) / 2.f;
    }
    //std::cout << "Found split: " << split << " (dim=" << dim << ") in " << link->core().min[dim] << " - " << link->core().max[dim] << std::endl;

    // subset and enqueue
    std::vector< std::vector<Point> > out_points(srp.out_link().size());
    for (size_t i = 0; i < (b->*points_).size(); ++i)
    {
      float x = (b->*points_)[i][dim];
      int loc = x < split ? 0 : 1;
      out_points[loc].push_back((b->*points_)[i]);
    }
    int pos = -1;
    for (int i = 0; i < k; ++i)
    {
      if (srp.out_link().target(i).gid == srp.gid())
      {
        (b->*points_).swap(out_points[i]);
        pos = i;
      }
      else
        srp.enqueue(srp.out_link().target(i), out_points[i]);
    }
//    float core_[2] = {link->core().min[dim], link->core().max[dim]}, 
//          regularcore_[2] = {link->regular_core().min[dim], link->regular_core().max[dim]},
//          bounds_[2] = {link->bounds().min[dim], link->bounds().max[dim]};
//    float prev_divs[2] = {link->divs().min[dim], link->divs().max[dim]};
    if (pos == 0) {
        link->divs().max[dim] = mid_div - 1;
        //link->regular_core().max[dim] = middle;
        link->core().max[dim] = split;
    }
    else {
        link->divs().min[dim] = mid_div;
        //link->regular_core().min[dim] = middle;
        link->core().min[dim] = split;
    }

/*    if (link->divs().min[dim] > link->divs().max[dim]) {
        fprintf(stderr, "dim=%d, prev_divs=[%f,%f], divs=[%f,%f]\n",
            dim, prev_divs[0], prev_divs[1], link->divs().min[dim], link->divs().max[dim]);
        assert (false);
    }
*/
/*    
    if (link->core().min[dim] == link->bounds().min[dim] && link->core().max[dim] == link->bounds().max[dim]) 
    {
      fprintf(stderr, "dim = %d, regularcore = [%f, %f] -> [%f, %f], core = [%f, %f] -> [%f, %f], bounds = [%f, %f] -> [%f, %f], split = %f\n", 
          dim, regularcore_[0], regularcore_[1], link->regular_core().min[dim], link->regular_core().max[dim], 
          core_[0], core_[1], link->core().min[dim], link->core().max[dim],
          bounds_[0], bounds_[1], link->bounds().min[dim], link->bounds().max[dim], split);
    }
*/
}

template<class Block, class Point>
void
diy::detail::ConstrainedKDTreePartition<Block,Point>::
dequeue_exchange(Block* b, const diy::ReduceProxy& srp, int dim) const
{
    int         lid  = srp.master()->lid(srp.gid());
    RCLink*     link = static_cast<RCLink*>(srp.master()->link(lid));

    for (int i = 0; i < srp.in_link().size(); ++i)
    {
      int nbr_gid = srp.in_link().target(i).gid;
      if (nbr_gid == srp.gid())
          continue;

      std::vector<Point>   in_points;
      srp.dequeue(nbr_gid, in_points);
      for (size_t j = 0; j < in_points.size(); ++j)
      {
        if (in_points[j][dim] < link->core().min[dim] || in_points[j][dim] > link->core().max[dim])
        {
            fprintf(stderr, "Warning: dequeued %f outside [%f,%f] (%d)\n",
                            in_points[j][dim], link->core().min[dim], link->core().max[dim], dim);
            std::abort();
        }
        (b->*points_).push_back(in_points[j]);
      }
    }
}

template<class Block, class Point>
void
diy::detail::ConstrainedKDTreePartition<Block,Point>::
update_neighbor_bounds(Bounds& bounds, float split, int dim, bool lower) const
{
    if (lower)
        bounds.max[dim] = split;
    else
        bounds.min[dim] = split;
}

template<class Block, class Point>
bool
diy::detail::ConstrainedKDTreePartition<Block,Point>::
intersects(const Bounds& x, const Bounds& y, int dim, const Bounds& domain) const
{
  return x.min[dim] <= y.max[dim] && y.min[dim] <= x.max[dim];
}

template<class Block, class Point>
float
diy::detail::ConstrainedKDTreePartition<Block,Point>::
find_split(const Bounds& changed, const Bounds& original) const // find_split(link->core(), link->bounds())
{
    for (int i = 0; i < dim_; ++i)
    {
        if (changed.min[i] != original.min[i])
            return changed.min[i];
        if (changed.max[i] != original.max[i])
            return changed.max[i];
    }

    assert(0);
    return -1;
}

#endif
