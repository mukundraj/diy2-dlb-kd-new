set (sources
  main.cpp
)

add_executable (convert ${sources})
target_link_libraries (convert
  ${NetCDF_LIBRARY}
  ${MPI_LIBRARIES}
)
