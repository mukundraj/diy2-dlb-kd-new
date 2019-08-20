set (sources
  tracer.cpp
  main.cpp
)

add_executable (sync_nek ${sources})
target_link_libraries (sync_nek
  appcommon dlb
)
