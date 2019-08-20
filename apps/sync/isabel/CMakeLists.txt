set (sources
  tracer.cpp
  main.cpp
)

add_executable (sync_isabel ${sources})
target_link_libraries (sync_isabel
  appcommon dlb
)
