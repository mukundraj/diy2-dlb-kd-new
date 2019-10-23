set (sources
  dlb.cc
  block.cpp
  cons_kdtree.cpp
  app.cpp
  app_sync.cpp
  misc.cpp
  #app_threadpool.cpp
)

set (proto
  config.proto
)

protobuf_generate_cpp (proto_sources proto_headers ${proto})

add_library (dlb ${sources} ${proto_sources} ${proto_headers})
target_link_libraries (dlb 
  async bil partitioner appcommon tinyxml
  ${PNETCDF_LIBRARY}
  ${MPI_LIBRARIES} 
  ${PROTOBUF_LIBRARIES}
  ${Boost_LIBRARIES}
)
