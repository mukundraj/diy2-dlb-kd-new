#include <cmath>
#include <algorithm>
#include <sstream>
#include <boost/foreach.hpp>
#include "tracer.h"
#include "common/advect.h"
#include "common/lerp.h"

const int max_trace_size = 2048;
const float stepsize = 0.2; // 450.f 
const int NUM_STEPS = 20;
const bool wrap = true;

CSyncGEOS5App::CSyncGEOS5App()
{
}

CSyncGEOS5App::~CSyncGEOS5App()
{
}

void CSyncGEOS5App::init(int argc, char **argv) {
  CPTApp_Sync::init(argc, argv);
}

void CSyncGEOS5App::initialize_particles(Block& b,
    std::vector<Particle>& particles) 
{
  // local-range analysis
  const int starts[4] = {59, 119, 0, 0},
            ends[4] = {96, 142, 71, 23};
  const int stride[] = {1, 1, 1, 1};
     
  for (int i = starts[0]; i <= ends[0]; i += stride[0]) 
    for (int j = starts[1]; j <= ends[1]; j += stride[1])
      for (int k = starts[2]; k <= ends[2]; k += stride[2]) 
        for (int l = starts[3]; l <= ends[3]; l += stride[3]) {
          const float idx[4] = {(float)i*1.f, (float)j*1.f, (float)k*1.f, (float)l*1.f}; 
          if (is_ptinblock(b, idx)) {
            Particle p;
            p.home_gid = b.gid;
            p[0] = i;  p[1] = j;  p[2] = k;  p[3] = l;
            particles.push_back(p);
          }  
        }
}
  
void CSyncGEOS5App::trace_particles(Block& b, 
    std::vector<Particle>& particles, 
    std::map<int, std::vector<Particle> >& unfinished_particles, 
    std::map<int, std::vector<Particle> >& finished_particles) 
{
  const float **vars = (const float**)(b.vars.data());
  int gst[4], gsz[4], lst[4], lsz[4];
  b.get_ghost_load_st_sz(num_dims(), gst, gsz, lst, lsz);

  BOOST_FOREACH (Particle& p, particles) {
    int steps = NUM_STEPS;
    while (p.num_steps < max_trace_size) {
      int rtn = trace_4D_geos5_rk1(gst, gsz, lst, lsz, vars, p.coords, wrap, stepsize);
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
      const int dst_gid = bound_gid(pt2gid(p.coords)); 
      unfinished_particles[dst_gid].push_back(p);
    }
  }
}

void CSyncGEOS5App::trace_particles_kdtree(Block& b, 
    std::vector<Particle>& particles, 
    std::map<int, std::vector<Particle> >& unfinished_particles, 
    std::map<int, std::vector<Particle> >& finished_particles) 
{
  const float **vars = (const float**)(b.vars.data());
  int gst[4], gsz[4], lst[4], lsz[4];
  b.get_ghost_load_st_sz(num_dims(), gst, gsz, lst, lsz);

  BOOST_FOREACH (Particle& p, particles) {
    int steps = NUM_STEPS;
    while (p.num_steps < max_trace_size) {
      int rtn = trace_4D_geos5_rk1(gst, gsz, lst, lsz, vars, p.coords, wrap, stepsize);
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

void CSyncGEOS5App::trace_particles_video(Block& b, 
      std::vector<Particle>& particles, 
      std::vector<Particle>& particles_first,
      std::vector<Particle>& particles_second
    )
{
  // TODO
}