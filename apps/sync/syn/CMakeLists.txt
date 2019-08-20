set (sources
  tracer.cpp
  main.cpp
)

add_executable (sync_syn ${sources})
target_link_libraries (sync_syn
  appcommon dlb
)