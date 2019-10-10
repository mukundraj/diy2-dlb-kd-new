#include "tracer.h"

using namespace std;

int main(int argc, char **argv)
{
  // MPI_Init(&argc, &argv);

  diy::mpi::environment env(argc, argv); // equivalent of MPI_Init(argc, argv)/MPI_Finalize()

  CSyncNekApp app;

  app.init(argc, argv);
  app.exec(); 
  app.deinit();

  // MPI_Finalize(); 
  return 0;
}