#include "tracer.h"

using namespace std;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  CSyncSynApp app;

  app.init(argc, argv);
  app.exec(); 
  app.deinit();

  MPI_Finalize(); 
  return 0;
}