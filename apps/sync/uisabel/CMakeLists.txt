set (sources
  tracer.cpp
  main.cpp
)

add_executable (sync_uisabel ${sources})
target_link_libraries (sync_uisabel
  appcommon dlb
)
