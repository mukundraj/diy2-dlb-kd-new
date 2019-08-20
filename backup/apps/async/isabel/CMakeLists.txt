set (sources
  tracer.cpp
  main.cpp
)

add_executable (async_isabel ${sources})
target_link_libraries (async_isabel
  appcommon dlb
)
