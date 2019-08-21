#include <cmath>
#include <algorithm>
#include <sstream>
#include <boost/foreach.hpp>
#include "tracer.h"
#include "common/advect.h"
#include "common/lerp.h"

const int max_trace_size = 2048;
const float stepsize = 1.0; 
const int NUM_STEPS = 200; // 50 for small data, 200 for 12GB data

CSyncNekApp::CSyncNekApp()
{
}

CSyncNekApp::~CSyncNekApp()
{
}

void CSyncNekApp::init(int argc, char **argv) {
  CPTApp_Sync::init(argc, argv);
}

void CSyncNekApp::initialize_particles(Block& b,
    std::vector<Particle>& particles) 
{
#if 1 // for k-d tree paper
    // {64, 128, 128} -> {256, 256, 256}
  const int stride[3] = {128, 128, 128};
  // const int stride[3] = {512, 512, 512};
  //const int stride[3] = {8, 8, 8};
  const float gap[3] = {1.f/(stride[0]-1) * (float)(domain_size()[0]-1),
                        1.f/(stride[1]-1) * (float)(domain_size()[1]-1),
                        1.f/(stride[2]-1) * (float)(domain_size()[2]-1)};

  float i = 0, j = 0, k = 0;
  while (i <= domain_size()[0]-1) {
    j =0; k = 0;
    while (j <= domain_size()[1]-1) {
      k = 0;
      while (k <= domain_size()[2]-1) {
        const float idx[3] = {i, j, k};
        if (is_ptinblock(b, idx)) {
          Particle p;
          p.home_gid = b.gid;
          p[0] = i;  p[1] = j;  p[2] = k;
          particles.push_back(p);
        }
        k += gap[2];
      }
      j += gap[1];
    }
    i += gap[0];
  }
#endif

#if 0 // original
  const int stride[3] = {4, 4, 4};
  for (int i = b.lb[0]; i <= b.ub[0]; i += stride[0])
    for (int j = b.lb[1]; j <= b.ub[1]; j += stride[1])
      for (int k = b.lb[2]; k <= b.ub[2]; k += stride[2]) {
        const int idx[3] = {i, j, k}; 
        const int id = idx2id(idx);
        Particle p;
        p.home_gid = b.gid;
        p[0] = i;  p[1] = j;  p[2] = k;
        particles.push_back(p);
      }
#endif
}
  
void CSyncNekApp::trace_particles(Block& b, 
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
      int rtn = trace_3D_rk1(gst, gsz, lst, lsz, vars, p.coords, stepsize);
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

#if 0
void CSyncNekApp::trace_particles_kdtree(Block& b, 
    std::vector<Particle>& particles, 
    std::map<int, std::vector<Particle> >& unfinished_particles_tbr, 
    std::map<int, std::vector<Particle> >& finished_particles)
{
  const float **vars = (const float**)(b.vars.data());
  //int gst[4], gsz[4];
  //b.get_ghost_st_sz(num_dims(), gst, gsz);
  int gst[4], gsz[4], lst[4], lsz[4];
  b.get_ghost_load_st_sz(num_dims(), gst, gsz, lst, lsz);

  std::vector<Particle> unfinished_particles;
  BOOST_FOREACH (Particle& p, b.particles) {
    int steps = NUM_STEPS;
    while (p.num_steps < max_trace_size) {
      //int rtn = trace_3D_rk1(gst, gsz, vars, p.coords, stepsize);
      int rtn = trace_3D_rk1(gst, gsz, lst, lsz, vars, p.coords, stepsize);
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
//fprintf(stderr, "p.num_steps = %d\n", p.num_steps);
    if (!inside_domain(p.coords) || p.num_steps >= max_trace_size)
      p.finished = true;

    if (p.finished) {
      //finished_particles[p.home_gid].push_back(p); // TODO
      b.num_particles_finished ++;
    } else {
      unfinished_particles.push_back(p);
    }
  }

  //b.num_particles = unfinished_particles.size();
  b.particles.clear();
  b.particles.insert(b.particles.end(), unfinished_particles.begin(), unfinished_particles.end());
}
#else
void CSyncNekApp::trace_particles_kdtree(Block& b, 
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
      int rtn = trace_3D_rk1(gst, gsz, lst, lsz, vars, p.coords, stepsize);
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
#endif

void CSyncNekApp::trace_particles_video(Block& b, 
      std::vector<Particle>& particles, 
      std::vector<Particle>& particles_first,
      std::vector<Particle>& particles_second
    )
{
  // TODO
}