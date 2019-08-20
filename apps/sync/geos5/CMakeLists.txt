set (sources
  tracer.cpp
  main.cpp
)

add_executable (sync_geos5 ${sources})
target_link_libraries (sync_geos5
  appcommon dlb
)
