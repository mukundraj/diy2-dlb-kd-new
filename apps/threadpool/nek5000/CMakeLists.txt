set (sources
  tracer.cpp
  main.cpp
)

add_executable (tp_nek5000 ${sources})
target_link_libraries (tp_nek5000
  appcommon dlb
)