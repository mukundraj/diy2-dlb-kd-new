#ifndef _SYNC_NEK_TRACER_H
#define _SYNC_NEK_TRACER_H

#include <algorithm>
#include <string>
#include <vector>
#include "src/app_sync.h"

class CSyncNekApp : public CPTApp_Sync {
public:
  CSyncNekApp();
  ~CSyncNekApp();

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

  void trace_particles_kdtree_predict(Block& b, int factor);

  void trace_particles_video(Block& b, 
      std::vector<Particle>&, 
      std::vector<Particle>&,
      std::vector<Particle>&
    );
  
};
#endif
