#ifndef _THREADPOOL_NEK_TRACER_H
#define _THREADPOOL_NEK_TRACER_H

#include <algorithm>
#include <string>
#include <vector>
#include "src/app_threadpool.h"

class CThreadPoolNekApp : public CPTApp_ThreadPool {
public:
  CThreadPoolNekApp();
  ~CThreadPoolNekApp();

  void init(int argc, char **argv);

public: 
  void initialize_particles(Block&,
      std::vector<Particle>&
    );
  
  void trace_particles(Block&,
      std::vector<Particle>&, // particles to trace
      std::map<int, std::vector<Particle> >&, // unfinished particles
      std::map<int, std::vector<Particle> >&  // finished particles
    );

  void trace_particles_kdtree(Block&,
      std::vector<Particle>&, // particles to trace
      std::map<int, std::vector<Particle> >&, // unfinished particles
      std::map<int, std::vector<Particle> >&  // finished particles
    );
  
};
#endif
