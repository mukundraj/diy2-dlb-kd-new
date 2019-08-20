#ifndef _SYNC_ISABEL_TRACER_H
#define _SYNC_ISABEL_TRACER_H

#include <algorithm>
#include <string>
#include <vector>
#include "src/app_sync.h"

class CSyncIsabelApp : public CPTApp_Sync {
public:
  CSyncIsabelApp();
  ~CSyncIsabelApp();

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
