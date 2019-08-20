#include <vector>
#include <cstdio>
#include <boost/foreach.hpp>
#include <diy/algorithms.hpp>
#include "dlb.h"
#include "cons_kdtree.h"
#include "partitioner/partitioner.h"

// constrained k-d tree decomposition with partial data duplication for SciVis'17 submission
// added by Jiang Zhang
// start from 12/06/2016

void populate_cons_kdtree_block(Block* blk, const diy::Master::ProxyWithLink& cp, void* aux_)
{ 
  PopulateMaster*        aux = (PopulateMaster*) aux_;
  //diy::Master* kdtree_master = (diy::Master*) aux;
  diy::Master* kdtree_master = aux->master;
  int                    dim = aux->dim;
  bool                   space_only = aux->space_only;

  CBounds domain;

  for (int i = 0; i < dim; i ++) {
    domain.min[i] = (blk->data_bounds).min[i];
    domain.max[i] = (blk->data_bounds).max[i];
  }

  ConstrainedKDTreeBlock* b = new ConstrainedKDTreeBlock;
  RCLink* l = new RCLink(dim, domain, domain);
  //l->set_regular_core(domain);

  CBounds divs;
  for (int i = 0; i < dim; i ++) {
    divs.min[i] = 0;
    divs.max[i] = aux->divs[i]-1;
  }
  l->set_divs(divs);

  kdtree_master->add(cp.gid(), b, l); // add a ConstrainedKDTreeBlock 
//fprintf(stderr, "before: dim = %d, ps= %d, divs=[%f,%f,%f], domain=[%f,%f,%f]\n", 
//  dim, blk->particles.size(), divs.max[0], divs.max[1], divs.max[2], 
//  domain.max[0], domain.max[1], domain.max[2]);
  // TODO determine which block for which b 
  // copy the particles over

  b->points.resize(blk->particles.size());
  for (size_t i = 0; i < blk->particles.size(); i ++) // TODO refine this
  {
    b->points[i][0] = blk->particles[i].coords[0];
    b->points[i][1] = blk->particles[i].coords[1];
    if (dim >= 3)
      b->points[i][2] = blk->particles[i].coords[2];
    if (space_only || dim == 4)
      b->points[i][3] = blk->particles[i].coords[3];
    b->points[i].home_gid = blk->particles[i].home_gid;
    //if (b->points[i].home_gid != blk->gid) assert (false);
    b->points[i].num_steps = blk->particles[i].num_steps;
    b->points[i].num_rounds = blk->particles[i].num_rounds;
    b->points[i].finished = blk->particles[i].finished;
  }

  //b->points.insert(b->points.end(), blk->particles.begin(), blk->particles.end());
  for (int i = 0; i < dim; i ++) {
    b->domain_mins[i] = (blk->data_bounds).min[i];
    b->domain_maxs[i] = (blk->data_bounds).max[i];
  }

//fprintf(stderr, "populate: dim = %d, blk->particles.size() = %d\n", dim, blk->particles.size());
}

void extract_cons_kdtree_block(ConstrainedKDTreeBlock* b, const diy::Master::ProxyWithLink& cp, void* aux_)
{
  ExtractMaster*      aux = (ExtractMaster*) aux_;
  diy::Master*  pt_master = aux->master;
  int                 dim = aux->dim;
  bool                space_only = aux->space_only;
  
  int pt_lid = pt_master->lid(cp.gid());

  Block* blk = (Block*) pt_master->block(pt_lid);    // assumes all the blocks are in memory

  // copy out the particles
  //blk->num_particles = b->points.size();
  blk->particles.clear();
  
  for (size_t i = 0; i < b->points.size(); i ++) { // TODO refine this
    Particle p;
    p.coords[0] = b->points[i][0];
    p.coords[1] = b->points[i][1];
    if (dim >= 3)
      p.coords[2] = b->points[i][2];
    if (space_only || dim == 4)
      p.coords[3] = b->points[i][3];
    //p.home_gid = b->points[i].home_gid;
    p.home_gid = blk->gid; // only for async parallel particle tracing (but sync one does not use this)
    p.num_steps = b->points[i].num_steps;
    p.num_rounds = b->points[i].num_rounds;
    p.finished = b->points[i].finished; //  = false
    blk->particles.push_back(p);
  }

  diy::RegularContinuousLink* kdtree_link = static_cast<diy::RegularContinuousLink*>(cp.link());

  for (int i = 0; i < dim; i ++) {
    blk->core_bounds.min[i] = kdtree_link->core().min[i];
    blk->core_bounds.max[i] = kdtree_link->core().max[i];
  }

  /*blk->particles.insert(blk->particles.end(), b->points.begin(), b->points.end());
  BOOST_FOREACH (Particle &p, blk->particles) {
    p.home_gid = blk->gid;
  }*/

  if (aux->async) {
    blk->after_kdtree = true;
    blk->num_particles_initialized = 0;
  }
//fprintf(stderr, "extract: dim = %d, b->points.size() = %d\n", dim, b->points.size());
  if (aux->first) {  // get lb, ub, glb, gub, lload, uload
    if (aux->constrained) {
      //diy::RegularContinuousLink* kdtree_link = static_cast<diy::RegularContinuousLink*>(cp.link());
      for (int i = 0; i < dim; i ++) {
        //blk->lb[i] = (int)(kdtree_link->regular_core().min[i]);
        //blk->ub[i] = (int)(kdtree_link->regular_core().max[i]) == b->domain_maxs[i] ? b->domain_maxs[i] : (int)(kdtree_link->regular_core().max[i]) - 1;
        blk->lb[i] = (int)(kdtree_link->divs().min[i]) * aux->block_size[i];
        blk->ub[i] = (int)(kdtree_link->divs().max[i]) == aux->divs[i]-1 ? b->domain_maxs[i] : (int)(kdtree_link->divs().max[i]+1)*aux->block_size[i]-1;

        blk->glb[i] = std::max(b->domain_mins[i], blk->lb[i] - aux->ghost_size[i]);
        blk->gub[i] = std::min(b->domain_maxs[i], blk->ub[i] + aux->ghost_size[i]);

        blk->lload[i] = blk->glb[i];
        blk->uload[i] = blk->gub[i];
      }
    }

    if (aux->async) {
      blk->after_kdtree = false;
      blk->num_particles_initialized = -1;
    }
  }

  delete b;     // safe to do since kdtree_master doesn't own the blocks (no create/destroy supplied)
}

double pt_cons_kdtree_exchange(
    diy::Master& master, const diy::Assigner& assigner, 
    std::vector<int> divisions, 
    const int dim, const bool space_only,
    int* block_size, int* ghost_size, 
    const bool constrained, const bool first, const bool async
  )
{
  diy::Master kdtree_master(master.communicator(), master.threads(), -1);
  //PopulateMaster populate_master = { &kdtree_master, dim, space_only };
  PopulateMaster populate_master;
  populate_master.master = &kdtree_master;
  populate_master.dim = dim;
  populate_master.space_only = space_only;
  for (int i = 0; i < dim; i ++)
    populate_master.divs[i] = divisions[i];
  master.foreach<Block>(&populate_cons_kdtree_block, &populate_master);
 
  CBounds domain;
  DBounds data_domain = master.block<Block>(master.loaded_block())->data_bounds;
  for (int i = 0; i < dim; i ++) {
    domain.min[i] = data_domain.min[i];
    domain.max[i] = data_domain.max[i];
  }
//fprintf(stderr, "start to k-d tree.\n");
  //if (sampling)
  //    diy::kdtree_sampling(kdtree_master, assigner, 3, domain, &KDTreeBlock::points, bins, wrap);
  size_t bins = 2048;      // 2048 histogram bins for uisabel and geos5; for nek small it is 8192
  diy::cons_kdtree(kdtree_master, assigner, dim, domain, divisions, block_size, ghost_size, constrained, first, &ConstrainedKDTreeBlock::points, bins);
//fprintf(stderr, "after k_d tree.\n");
  //diy::Master pt_master(master.communicator(),  master.threads());//, -1);
  //int gs[4] = {ghost_size[0], ghost_size[1], ghost_size[2], ghost_size[3]};
  ExtractMaster extract_master;// = { &pt_master, gs };
  extract_master.master = &master;
  extract_master.dim = dim;
  extract_master.space_only = space_only;
  for (int i = 0; i < dim; i ++) {
    extract_master.divs[i] = divisions[i];
    extract_master.block_size[i] = block_size[i];
    extract_master.ghost_size[i] = ghost_size[i];
  }
  extract_master.constrained = constrained;
  extract_master.first = first;
  extract_master.async = async;
  kdtree_master.foreach<ConstrainedKDTreeBlock>(&extract_cons_kdtree_block, &extract_master);
  //kdtree_master.foreach<ConstrainedKDTreeBlock>(&extract_cons_kdtree_block, &master);
  master.set_expected(kdtree_master.expected());

  return 0;
}
