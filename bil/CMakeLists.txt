set (BIL_SOURCES 
  bil.c
  bil_block.c
  bil_misc.c
  bil_pio.c
  bil_sched.c
  bil_timing.c
)

add_library (bil ${BIL_SOURCES})
target_link_libraries (bil ${MPI_LIBRARIES} ${PNETCDF_LIBRARY})
