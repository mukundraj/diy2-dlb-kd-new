#include <cmath>
#include <algorithm>
#include <sstream>
#include <boost/foreach.hpp>
#include "tracer.h"
#include "common/advect.h"
#include "common/lerp.h"

const int max_trace_size = 1024;
const float stepsize = 450.f; // 450.f in k-d tree paper
const int NUM_STEPS = 50; // 50 in k-d tree paper

CSyncUIsabelApp::CSyncUIsabelApp()
{
}

CSyncUIsabelApp::~CSyncUIsabelApp()
{
}

void CSyncUIsabelApp::init(int argc, char **argv) {
  CPTApp_Sync::init(argc, argv);
}

void CSyncUIsabelApp::initialize_particles(Block& b,
    std::vector<Particle>& particles) 
{
#if 0 // for the comparison between dynamic data repartitioning and k-d tree
  const int stride[4] = {192, 192, 32, 48};
  const float gap[4] = {1.f/(stride[0]-1) * (float)(domain_size()[0]-1),
                        1.f/(stride[1]-1) * (float)(domain_size()[1]-1),
                        1.f/(stride[2]-1) * (float)(domain_size()[2]-1),
                        1.f/(stride[3]-1) * (float)(domain_size()[3]-1)};

  float i = 0, j = 0, k = 0, l = 0;
  while (i <= domain_size()[0]-1) {
    j = 0; k = 0; l = 0;
    while (j <= domain_size()[1]-1) {
      k = 0; l = 0;
      while (k <= domain_size()[2]-1) {
        l = 0;
        while (l <= domain_size()[3]-1) {
          const float idx[4] = {i, j, k, l};
          if (is_ptinblock(b, idx)) {
            Particle p;
            p.home_gid = b.gid;
            p[0] = i;  p[1] = j;  p[2] = k; p[3] = l;
            particles.push_back(p);
          }
          l += gap[3];
        }
        k += gap[2];
      }
      j += gap[1];
    }
    i += gap[0];
  }
#endif 

#if 1 // used in k-d tree paper
  // {16, 16, 12, 48} -> {64, 64, 24, 48}
  //const int stride[4] = {64, 64, 24, 48};
  const int stride[4] = {128, 128, 32, 48};
  const float gap[4] = {1.f/(stride[0]-1) * (float)(domain_size()[0]-1),
                        1.f/(stride[1]-1) * (float)(domain_size()[1]-1),
                        1.f/(stride[2]-1) * (float)(domain_size()[2]-1),
                        1.f/(stride[3]-1) * (float)(domain_size()[3]-1)};

  float i = 0, j = 0, k = 0, l = 0;
  while (i <= domain_size()[0]-1) {
    j = 0; k = 0; l = 0;
    while (j <= domain_size()[1]-1) {
      k = 0; l = 0;
      while (k <= domain_size()[2]-1) {
        l = 0;
        while (l <= domain_size()[3]-1) {
          const float idx[4] = {i, j, k, l};
          if (is_ptinblock(b, idx)) {
            Particle p;
            p.home_gid = b.gid;
            p[0] = i;  p[1] = j;  p[2] = k; p[3] = l;
            particles.push_back(p);
          }
          l += gap[3];
        }
        k += gap[2];
      }
      j += gap[1];
    }
    i += gap[0];
  }
#endif 

#if 0 // original
  const int stride[4] = {2, 2, 2, 1};
  for (int i = b.lb[0]; i <= b.ub[0]; i += stride[0])
    for (int j = b.lb[1]; j <= b.ub[1]; j += stride[1])
      for (int k = b.lb[2]; k <= b.ub[2]; k += stride[2]) 
        for (int l = b.lb[3]; l <= b.ub[3]; l += stride[3]) {
          const int idx[4] = {i, j, k, l}; 
          const int id = idx2id(idx);
          Particle p;
          p.id = id; 
          p.home_gid = b.gid;
          p[0] = i;  p[1] = j;  p[2] = k; p[3] = l;
          particles.push_back(p);
        }
#endif
}
  
void CSyncUIsabelApp::trace_particles(Block& b, 
    std::vector<Particle>& particles, 
    std::map<int, std::vector<Particle> >& unfinished_particles, 
    std::map<int, std::vector<Particle> >& finished_particles
  )  // TODO no need to make sure destination
{
  const float **vars = (const float**)(b.vars.data());
  int gst[4], gsz[4], lst[4], lsz[4];
  b.get_ghost_load_st_sz(num_dims(), gst, gsz, lst, lsz);

  BOOST_FOREACH (Particle& p, particles) {
    int steps = NUM_STEPS;
    while (p.num_steps < max_trace_size) {
      int rtn = trace_4D_isabel_rk1(gst, gsz, lst, lsz, vars, p.coords, stepsize);
      if (rtn == TRACE_OUT_OF_BOUND) break; // out of ghost size
      add_workload();
      if (rtn == TRACE_CRITICAL_POINT || rtn == TRACE_NO_VALUE) {
        p.finished = true;
        break;
      }
      p.num_steps ++;
      steps --;
      
      if (steps <= 0) break;
    }

    if (!inside_domain(p.coords) || p.num_steps >= max_trace_size)
      p.finished = true;

    if (p.finished) {
      finished_particles[p.home_gid].push_back(p);
      _local_done ++;
    } else {
      const int dst_gid = bound_gid(pt2gid(p.coords)); // TODO
      unfinished_particles[dst_gid].push_back(p);
    }
  }
}

void CSyncUIsabelApp::trace_particles_kdtree(Block& b,
    std::vector<Particle>& particles, 
    std::map<int, std::vector<Particle> >& unfinished_particles, 
    std::map<int, std::vector<Particle> >& finished_particles
  ) 
{
  const float **vars = (const float**)(b.vars.data());
  int gst[4], gsz[4], lst[4], lsz[4];
  b.get_ghost_load_st_sz(num_dims(), gst, gsz, lst, lsz);

  BOOST_FOREACH (Particle& p, particles) {
    int steps = NUM_STEPS;
    while (p.num_steps < max_trace_size) {
      int rtn = trace_4D_isabel_rk1(gst, gsz, lst, lsz, vars, p.coords, stepsize);
      if (rtn == TRACE_OUT_OF_BOUND) break; // out of ghost size

      add_workload();
      if (rtn == TRACE_CRITICAL_POINT || rtn == TRACE_NO_VALUE) {
        p.finished = true;
        break;
      }
      p.num_steps ++;
      steps --;
      
      if (steps <= 0) break;
    }

    if (!inside_domain(p.coords) || p.num_steps >= max_trace_size)
      p.finished = true;

    if (p.finished) {
      _local_done ++;
    } else {
      b.particles.push_back(p);
    }
  }
}


void CSyncUIsabelApp::trace_particles_video(Block& b, 
      std::vector<Particle>& particles, 
      std::vector<Particle>& particles_first,
      std::vector<Particle>& particles_second
    )
{
  // TODO
}