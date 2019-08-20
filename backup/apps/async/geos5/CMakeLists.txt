set (sources
  tracer.cpp
  main.cpp
)

add_executable (async_geos5 ${sources})
target_link_libraries (async_geos5
  appcommon dlb
)
