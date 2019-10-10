#ifndef _DLB_H
#define _DLB_H

#include <cstdio>
#include <vector>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cmath>
#include <climits>
#include <stdint.h>
#include <diy/serialization.hpp>

// NetCDF error handling
#define NC_SAFE_CALL(call) {\
  int retval = call;\
  if (retval != 0) {\
    fprintf(stderr, "[NetCDF Error] %s, in file '%s', line %i.\n", nc_strerror(retval), __FILE__, __LINE__); \
    exit(EXIT_FAILURE); \
  }\
}

// Parallel-NetCDF error handling
#define PNC_SAFE_CALL(call) {\
  int retval = call;\
  if (retval != 0) {\
      fprintf(stderr, "[PNetCDF Error] %s, in file '%s', line %i.\n", ncmpi_strerror(retval), __FILE__, __LINE__); \
      exit(EXIT_FAILURE); \
  }\
}

struct Particle {
  int id; 
  int home_gid;
  int num_steps;
  int num_rounds;
  float coords[4];
  bool finished;
  int wgt;
  bool real;
  int num_esteps; // number of epoch steps
  bool epoch_finished;

  Particle() {
    id = 0;
    num_steps = 0;
    num_esteps = 0;
    num_rounds = 0;
    memset(coords, sizeof(float)*4, 0);
    finished = false;
    wgt = 1;
    real = true;
    epoch_finished = false;

  }

  const float& operator [](int i) const {return coords[i];}
  float& operator [](int i) {return coords[i];}
};

// namespace diy {
//   template <> struct Serialization<Particle> {
//     static void save(BinaryBuffer& bb, const Particle& m) {
//       diy::save(bb, m.id);
//       diy::save(bb, m.home_gid);
//       diy::save(bb, m.num_steps);
//       diy::save(bb, m.num_rounds);
//       diy::save(bb, m.coords, 4);
//       diy::save(bb, m.finished);
//       diy::save(bb, m.wgt);
//       diy::save(bb, m.real);
//       diy::save(bb, m.num_esteps);
//     }
    
//     static void load(BinaryBuffer& bb, Particle& m) {
//       diy::load(bb, m.id);
//       diy::load(bb, m.home_gid);
//       diy::load(bb, m.num_steps);
//       diy::load(bb, m.num_rounds);
//       diy::load(bb, m.coords, 4);
//       diy::load(bb, m.finished);
//       diy::load(bb, m.wgt);
//       diy::load(bb, m.real);
//       diy::load(bb, m.num_esteps);
//     }
//   };
// }


struct Message {
  enum {
    MESSAGE_TRACE = 0,
    MESSAGE_INIT,
    MESSAGE_FINISH
  };

  unsigned char type;
  int src_gid, dst_gid; // mandatory. TODO: change to unsigned short (max=65536)
  double src_workload;
  uint64_t num_finished_particles;

  std::vector<Particle> unfinished_particles;
  std::vector<Particle> finished_particles;

  int size() const {return unfinished_particles.size() + finished_particles.size();}
  void swap(Message& m) {
    type = m.type;
    src_gid = m.src_gid; 
    dst_gid = m.dst_gid;
    src_workload = m.src_workload;
    num_finished_particles = m.num_finished_particles;
    unfinished_particles.swap(m.unfinished_particles);
    finished_particles.swap(m.finished_particles);
  }

  Message() : 
    type(MESSAGE_TRACE),
    src_gid(0), dst_gid(0), src_workload(0), num_finished_particles(0)
  {}
};

namespace diy {
  template <> struct Serialization<Message> {
    static void save(BinaryBuffer& bb, const Message& m) {
      diy::save(bb, m.src_gid);
      diy::save(bb, m.dst_gid);
      diy::save(bb, m.type);
      diy::save(bb, m.src_workload);
      diy::save(bb, m.num_finished_particles);
      diy::save(bb, m.unfinished_particles);
      diy::save(bb, m.finished_particles);
    }
    
    static void load(BinaryBuffer& bb, Message& m) {
      diy::load(bb, m.src_gid);
      diy::load(bb, m.dst_gid);
      diy::load(bb, m.type);
      diy::load(bb, m.src_workload);
      diy::load(bb, m.num_finished_particles);
      diy::load(bb, m.unfinished_particles);
      diy::load(bb, m.finished_particles);
    }
  };
}

/*
struct Message {
  int src_gid, dst_gid; // mandatory. TODO: change to unsigned short (max=65536)
  uint64_t num_finished_particles;
  std::vector<Particle> unfinished_particles;
  //std::vector<Particle> finished_particles;

  int size() const { return unfinished_particles.size(); }
  void swap(Message& m) {
    src_gid = m.src_gid; 
    dst_gid = m.dst_gid;
    num_finished_particles = m.num_finished_particles;
    unfinished_particles.swap(m.unfinished_particles);
    //finished_particles.swap(m.finished_particles);
  }

  Message() : 
    src_gid(0), dst_gid(0), num_finished_particles(0)
  {}
};

namespace diy {
  template <> struct Serialization<Message> {
    static void save(BinaryBuffer& bb, const Message& m) {
      diy::save(bb, m.src_gid);
      diy::save(bb, m.dst_gid);
      diy::save(bb, m.num_finished_particles);
      diy::save(bb, m.unfinished_particles);
      //diy::save(bb, m.finished_particles);
    }
    
    static void load(BinaryBuffer& bb, Message& m) {
      diy::load(bb, m.src_gid);
      diy::load(bb, m.dst_gid);
      diy::load(bb, m.num_finished_particles);
      diy::load(bb, m.unfinished_particles);
      //diy::load(bb, m.finished_particles);
    }
  };
}
*/

template <typename T>
static double imbalance(std::vector<T> loads)
{
  double l_max=0, l_ave=0;
  for (int i=0; i<loads.size(); i++) {
    l_max = std::max(l_max, (double)loads[i]);
    l_ave += loads[i];
  }
  l_ave = l_ave/loads.size();
  double lambda = (l_max - l_ave) / l_max;
  return lambda; 
}

std::string pt2buf(const int num_dims, const float *pt);
void buf2pt(const int num_dims, const std::string &buf, float *pt);

std::string id2str(const int num_dims, const int *id);
void str2id(const int num_dims, const std::string &buf, int *id);

void fill_bounds(
    const int& ghost_block_size, 
    const int* domain_minmax, 
    const int* core, 
    int* bounds
  );

#endif
