#include "tracer.h"

using namespace std;

int main(int argc, char **argv)
{
  //MPI_Init(&argc, &argv);
  const int required = CPTApp_ThreadPool::mpi_thread_level();
  int provided;

  MPI_Init_thread(&argc, &argv, required, &provided);
  assert(required == provided);

  CThreadPoolNekApp app;

  app.init(argc, argv);
  app.exec(); 
  app.deinit();

  MPI_Finalize(); 
  return 0;
}