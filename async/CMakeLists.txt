set (async_sources
  comm.cpp
  async_comm.cpp
  status.cpp
)

add_library (async ${async_sources})
target_link_libraries (async ${MPI_LIBRARIES})
