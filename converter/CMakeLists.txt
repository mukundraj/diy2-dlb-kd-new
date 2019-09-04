set (sources
  main.cpp
)

add_executable (convert ${sources})
target_link_libraries (convert
  ${NetCDF_LIBRARY}
  ${MPI_LIBRARIES}
)


add_executable (subset_selector subset_selector.cpp)
target_link_libraries (subset_selector
  ${NetCDF_LIBRARY}
  ${MPI_LIBRARIES}
)