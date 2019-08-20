#include "partitioner.h"
#include <cstdio>
#include <math.h>
#include <iostream>
#include <algorithm>

using namespace std;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  const int domain_size[3] = {500, 500, 100};
  const int  ghost_size[3] = {1, 1, 1};

  RegularPartitioner partitioner(comm, 3, 1, 3);
  partitioner.set_domain_size(domain_size);
  partitioner.set_num_blocks(size);
  partitioner.update_partition();

  partitioner.set_ghost(ghost_size);
  partitioner.update_ghosts(false);

  partitioner.print_blk_info();

  MPI_Finalize();
  return 0;
}
