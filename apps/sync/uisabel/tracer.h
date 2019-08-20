#ifndef _SYNC_UISABEL_TRACER_H
#define _SYNC_UISABEL_TRACER_H

#include <algorithm>
#include <string>
#include <vector>
#include "src/app_sync.h"

class CSyncUIsabelApp : public CPTApp_Sync {
public:
  CSyncUIsabelApp();
  ~CSyncUIsabelApp();

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

  void trace_particles_video(Block& b, 
      std::vector<Particle>&, 
      std::vector<Particle>&,
      std::vector<Particle>&
    );

};

#endif
