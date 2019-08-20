#include <cmath>
#include <algorithm>
#include <sstream>
#include <boost/foreach.hpp>
#include "tracer.h"
#include "common/advect.h"
#include "common/lerp.h"

const int max_trace_size = 2000;
const float stepsize = 0.2; // 450.f 
const bool wrap = true;

CAsyncGEOS5App::CAsyncGEOS5App()
{
}

CAsyncGEOS5App::~CAsyncGEOS5App()
{
}

void CAsyncGEOS5App::init(int argc, char **argv) {
  CPTApp_Async::init(argc, argv);
}

#if 1
void CAsyncGEOS5App::initialize_particles(Block& b,
    std::vector<Particle>& particles) 
{
#if 0 // full-range
  // {16, 9, 9, 24} -> {64, 36, 18, 24}
  //const int stride[4] = {64, 36, 18, 24};
  const int stride[4] = {64, 36, 18, 24};
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
#else  // local range
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
/*
  const int starts[4] = {59, 119, 0, 0},
            ends[4] = {96, 142, 71, 23};
  const int stride[] = {38, 24, 72, 24};

  const float gap[4] = {1.f/(stride[0]-1) * (float)(ends[0] - starts[0]),
                        1.f/(stride[1]-1) * (float)(ends[1] - starts[1]),
                        1.f/(stride[2]-1) * (float)(ends[2] - starts[2]),
                        1.f/(stride[3]-1) * (float)(ends[3] - starts[3])};

  float i = starts[0], j = starts[1], k = starts[2], l = starts[3];
  while (i <= ends[0]) {
    j = starts[1], k = starts[2], l = starts[3];
    while (j <= ends[1]) {
      k = starts[2], l = starts[3];
      while (k <= ends[2]) {
        l = starts[3];
        while (l <= ends[3]) {
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
*/
#endif
} 
#else
void CAsyncGEOS5App::initialize_particles(Block& b,
    std::vector<Particle>& particles) 
{
#if 1 // full-range analysis
  const int stride[4] = {8, 4, 4, 1};
  for (int i = 0; i <= domain_size()[0]; i += stride[0])
    for (int j = 0; j <= domain_size()[1]; j += stride[1])
      for (int k = 0; k <= domain_size()[2]; k += stride[2])
        for (int l = 0; l <= domain_size()[3]; l += stride[3]) {
          const float idx[4] = {(float)i*1.f, (float)j*1.f, (float)k*1.f, (float)l*1.f}; 
          if (is_ptinblock(b, idx)) {
            Particle p;
            p.home_gid = b.gid;
            p[0] = i;  p[1] = j;  p[2] = k; p[3] = l;
            particles.push_back(p);
          }
        }
#else // local-range analysis 
  const int starts[3] = {59, 119, 0},
            ends[3] = {96, 142, 71};
  const int stride[] = {2, 2, 2};

  /*
  const int starts[3] = {188, 116, 0},
            ends[3] = {234, 156, 40};
  const int stride[] = {2, 2, 4};
  */
     
  for (int i = starts[0]; i < ends[0]; i += stride[0]) 
    for (int j = starts[1]; j < ends[1]; j += stride[1])
      for (int k = starts[2]; k < ends[2]; k += stride[2]) {
        const int idx[3] = {i, j, k};
        if (inside_lb_ub(3, b.lb, b.ub, idx)) {
          const int id = idx2id(idx);
          Particle p;
          p.id = id; 
          p.home_gid = b.gid;
          p[0] = i;  p[1] = j;  p[2] = k;
          particles.push_back(p);
        }
      }
#endif
}
#endif
/*
void CAsyncGEOS5App::initialize_particles(Block& b,
    std::vector<Particle>& particles) 
{
#if 1 // full-range analysis
  const int stride[4] = {8, 4, 4, 1};
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
#else // local-range analysis 
  const int starts[3] = {59, 119, 0},
            ends[3] = {96, 142, 71};
  const int stride[] = {2, 2, 2};
     
  for (int i = starts[0]; i < ends[0]; i += stride[0]) 
    for (int j = starts[1]; j < ends[1]; j += stride[1])
      for (int k = starts[2]; k < ends[2]; k += stride[2]) {
        const int idx[3] = {i, j, k};
        if (inside_lb_ub(3, b.lb, b.ub, idx)) {
          const int id = idx2id(idx);
          Particle p;
          p.id = id; 
          p.home_gid = b.gid;
          p[0] = i;  p[1] = j;  p[2] = k;
          particles.push_back(p);
        }
      }
#endif
}
*/  
void CAsyncGEOS5App::trace_particles(Block& b, 
    std::vector<Particle>& particles, 
    std::map<int, std::vector<Particle> >& unfinished_particles, 
    std::map<int, std::vector<Particle> >& finished_particles
  ) 
{
  const float **vars = (const float**)(b.vars.data());
  int gst[4], gsz[4], lst[4], lsz[4];
  b.get_ghost_load_st_sz(num_dims(), gst, gsz, lst, lsz);

  BOOST_FOREACH (Particle& p, particles) {
    while (p.num_steps < max_trace_size) {
      int rtn = trace_4D_geos5_rk1(gst, gsz, lst, lsz, vars, p.coords, wrap, stepsize);
      if (rtn == TRACE_OUT_OF_BOUND) break; // out of ghost size
      add_workload();
      if (rtn == TRACE_CRITICAL_POINT || rtn == TRACE_NO_VALUE) {
        p.finished = true;
        break;
      }
      p.num_steps ++;
    }

    if (!inside_domain(p.coords) || p.num_steps >= max_trace_size)
      p.finished = true;

    if (p.finished) {
      finished_particles[p.home_gid].push_back(p);
    } else {
      const int dst_gid = pt2gid(p.coords); // TODO
      unfinished_particles[dst_gid].push_back(p);
    }
  }
}
/*
void CAsyncGEOS5App::trace_particles_kdtree(Block& b, 
    std::vector<Particle>& particles, 
    std::map<int, std::vector<Particle> >& unfinished_particles, 
    std::map<int, std::vector<Particle> >& finished_particles
  ) 
{
  const float **vars = (const float**)(b.vars.data());
  int gst[4], gsz[4], lst[4], lsz[4];
  b.get_ghost_load_st_sz(num_dims(), gst, gsz, lst, lsz);

  int round_trace_size = max_trace_size / NUM_ROUNDS;
  BOOST_FOREACH (Particle& p, particles) {
    int num_steps = p.num_rounds * round_trace_size;
    while (p.num_steps - num_steps < round_trace_size) {
      int rtn = trace_4D_geos5_rk1(gst, gsz, lst, lsz, vars, p.coords, wrap, stepsize);
      if (rtn == TRACE_OUT_OF_BOUND) break; // out of ghost size
      add_workload();
      if (rtn == TRACE_CRITICAL_POINT || rtn == TRACE_NO_VALUE) {
        p.finished = true;
        break;
      }
      p.num_steps ++;
    }

    if (!inside_domain(p.coords) || p.num_steps >= max_trace_size)
      p.finished = true;

    if (p.finished) {
      finished_particles[p.home_gid].push_back(p);
    } else {
      if (p.num_steps - num_steps == round_trace_size) { // round finished
        finished_particles[p.home_gid].push_back(p); // not really finished
        
        p.home_gid = b.gid;
        p.num_rounds ++;
        b.particles.push_back(p);
      } else { // reach the boundary, to be sent to other processes
        const int dst_gid = bound_gid(pt2gid(p.coords)); // TODO
        unfinished_particles[dst_gid].push_back(p);
      }
    }
  }
}
*/

void CAsyncGEOS5App::trace_particles_kdtree(Block& b, 
    std::vector<Particle>& particles, 
    std::map<int, std::vector<Particle> >& unfinished_particles, 
    std::map<int, int>& num_finished_particles
  ) 
{
  const float **vars = (const float**)(b.vars.data());
  int gst[4], gsz[4], lst[4], lsz[4];
  b.get_ghost_load_st_sz(num_dims(), gst, gsz, lst, lsz);

  int round_trace_size = max_trace_size / NUM_ROUNDS;
  BOOST_FOREACH (Particle& p, particles) {
    int num_steps = p.num_rounds * round_trace_size;
    while (p.num_steps - num_steps < round_trace_size) {
      int rtn = trace_4D_geos5_rk1(gst, gsz, lst, lsz, vars, p.coords, wrap, stepsize);
      if (rtn == TRACE_OUT_OF_BOUND) break; // out of ghost size
      add_workload();
      if (rtn == TRACE_CRITICAL_POINT || rtn == TRACE_NO_VALUE) {
        p.finished = true;
        break;
      }
      p.num_steps ++;
    }

    if (!inside_domain(p.coords) || p.num_steps >= max_trace_size)
      p.finished = true;

    if (p.finished) {
      num_finished_particles[p.home_gid] ++;
    } else {
      if (p.num_steps - num_steps == round_trace_size) { // round finished
        num_finished_particles[p.home_gid] ++; // not really finished
        
        p.home_gid = b.gid;
        p.num_rounds ++;
        b.particles.push_back(p); // b.kd_particles.push_back(p);
      } else { // reach the boundary, to be sent to other processes
        const int dst_gid = bound_gid(pt2gid(p.coords)); // TODO
        unfinished_particles[dst_gid].push_back(p);
      }
    }
  }
}