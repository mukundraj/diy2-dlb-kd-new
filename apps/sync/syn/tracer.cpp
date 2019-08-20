#include <cmath>
#include <algorithm>
#include <sstream>
#include <boost/foreach.hpp>
#include "tracer.h"
#include "common/advect.h"
#include "common/lerp.h"

const int max_trace_size = 1024;
const float stepsize = 0.2f; 
const int NUM_STEPS = 90; 

CSyncSynApp::CSyncSynApp()
{
}

CSyncSynApp::~CSyncSynApp()
{
}

void CSyncSynApp::init(int argc, char **argv) {
  CPTApp_Sync::init(argc, argv);
}

void CSyncSynApp::initialize_particles(Block& b,
    std::vector<Particle>& particles) 
{
#if 1
    // {64, 128, 128} -> {256, 256, 256}
  const int stride[2] = {40, 40};
  //const int stride[3] = {8, 8, 8};
  const float gap[2] = {1.f/(stride[0]-1) * (float)(domain_size()[0]-1),
                        1.f/(stride[1]-1) * (float)(domain_size()[1]-1)};

  float i = 0, j = 0;
  while (i <= domain_size()[0]-1) {
    j = 0; 
    while (j <= domain_size()[1]-1) {
      const float idx[2] = {i, j};
      if (is_ptinblock(b, idx) && i!=0 && j!=0) {
        Particle p;
        p.home_gid = b.gid;
        p[0] = i;  p[1] = j; 
        particles.push_back(p);
      }
      j += gap[1];
    }
    i += gap[0];
  }
#else
  const int stride[2] = {4, 4};
  for (int i = b.lb[0]; i <= b.ub[0]; i += stride[0])
    for (int j = b.lb[1]; j <= b.ub[1]; j += stride[1]) {
      const int idx[2] = {i, j}; 
      const int id = idx2id(idx);
      Particle p;
      p.home_gid = b.gid;
      p[0] = i;  p[1] = j;
      particles.push_back(p);
    }
#endif
}
  
void CSyncSynApp::trace_particles(Block& b, 
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
      int rtn = trace_2D_rk1(gst, gsz, lst, lsz, vars, p.coords, stepsize);
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

void CSyncSynApp::trace_particles_kdtree(Block& b, 
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
      int rtn = trace_2D_rk1(gst, gsz, lst, lsz, vars, p.coords, stepsize);
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

void CSyncSynApp::trace_particles_video(Block& b, 
    std::vector<Particle>& particles, 
    std::vector<Particle>& particles_first,
    std::vector<Particle>& particles_second)
{
  const float **vars = (const float**)(b.vars.data());
  int gst[4], gsz[4], lst[4], lsz[4];
  b.get_ghost_load_st_sz(num_dims(), gst, gsz, lst, lsz);

  BOOST_FOREACH (Particle& p, particles) {
    int steps = NUM_STEPS;
    while (p.num_steps < max_trace_size) {
      int rtn = trace_2D_rk1(gst, gsz, lst, lsz, vars, p.coords, stepsize);
      if (rtn == TRACE_OUT_OF_BOUND) break; // out of ghost size
      add_workload();
      if (rtn == TRACE_CRITICAL_POINT || rtn == TRACE_NO_VALUE) {
        p.finished = true;
        break;
      }
      p.num_steps ++;
      steps --;
      
      if (steps == 60) 
        particles_first.push_back(p);
      if (steps == 30) 
        particles_second.push_back(p);

      if (steps <= 0) break;
    }

    if (!inside_domain(p.coords) || p.num_steps >= max_trace_size)
      p.finished = true;

    if (p.finished) {
      _local_done ++;
    } else {
      b.particles.push_back(p);

      if (steps >= 60) {
        particles_first.push_back(p);
        particles_second.push_back(p);
      }
      if (steps < 60 && steps >= 30) 
        particles_second.push_back(p);
    }
  }
}

