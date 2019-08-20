set (sources
  partitioner.cpp
)

add_library (partitioner ${sources})
target_link_libraries (partitioner 
  bil
  ${MPI_LIBRARIES} 
)

set (test_sources
  testpartition.cpp
)

add_executable (testpartition ${test_sources})
target_link_libraries(testpartition partitioner)
