#ifndef _ASYNC_UISABEL_TRACER_H
#define _ASYNC_UISABEL_TRACER_H

#include <algorithm>
#include <string>
#include <vector>
#include "src/app_async.h"

class CAsyncUIsabelApp : public CPTApp_Async {
public:
  CAsyncUIsabelApp();
  ~CAsyncUIsabelApp();

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
/*
  void trace_particles_kdtree(Block&,
      std::vector<Particle>&, // particles to trace
      std::map<int, std::vector<Particle> >&, // unfinished particles
      std::map<int, std::vector<Particle> >&  // finished particles
    );
*/
  void trace_particles_kdtree(Block&,
      std::vector<Particle>&, // particles to trace
      std::map<int, std::vector<Particle> >&, // unfinished particles
      std::map<int, int>&  // number of finished particles
    );

};

#endif
