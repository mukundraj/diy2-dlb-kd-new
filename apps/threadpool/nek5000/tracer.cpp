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

CThreadPoolNekApp::CThreadPoolNekApp()
{
}

CThreadPoolNekApp::~CThreadPoolNekApp()
{
}

void CThreadPoolNekApp::init(int argc, char **argv) {
  CPTApp_ThreadPool::init(argc, argv);
}

void CThreadPoolNekApp::initialize_particles(Block& b,
    std::vector<Particle>& particles) 
{
  // {64, 128, 128} -> {256, 256, 256}
  const int stride[3] = {512, 512, 512};
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
}
  
void CThreadPoolNekApp::trace_particles(Block& b, 
    std::vector<Particle>& particles, 
    std::map<int, std::vector<Particle> >& unfinished_particles, 
    std::map<int, std::vector<Particle> >& finished_particles) 
{
  const float **vars = (const float**)(b.vars.data());
  int workload = 0;
  int gst[4], gsz[4], lst[4], lsz[4];
  b.get_ghost_load_st_sz(num_dims(), gst, gsz, lst, lsz);

  BOOST_FOREACH (Particle& p, particles) {
    int steps = NUM_STEPS;
    while (p.num_steps < max_trace_size) {
      int rtn = trace_3D_rk1(gst, gsz, lst, lsz, vars, p.coords, stepsize);
      if (rtn == TRACE_OUT_OF_BOUND) break; // out of ghost size
      workload ++;
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
    } else {
      const int dst_gid = bound_gid(pt2gid(p.coords)); // TODO
      unfinished_particles[dst_gid].push_back(p);
    }
  }

  this->inc_workload(b, workload);
}

void CThreadPoolNekApp::trace_particles_kdtree(Block& b, 
    std::vector<Particle>& particles, 
    std::map<int, std::vector<Particle> >& unfinished_particles, 
    std::map<int, std::vector<Particle> >& finished_particles)
{
  // TODO
}
