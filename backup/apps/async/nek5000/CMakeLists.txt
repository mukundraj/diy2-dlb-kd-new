set (sources
  tracer.cpp
  main.cpp
)

add_executable (async_nek ${sources})
target_link_libraries (async_nek
  appcommon dlb
)
