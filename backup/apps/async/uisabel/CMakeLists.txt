set (sources
  tracer.cpp
  main.cpp
)

add_executable (async_uisabel ${sources})
target_link_libraries (async_uisabel
  appcommon dlb
)
