#undef NDEBUG
#include <iostream>
#include <vector>
#include <mpi.h>
#include <cassert>

#include "util/subcomm.h"

#define EV(x) (#x) << "=" << (x) << " "

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  MPI_Comm comm_world;
  MPI_Comm comm_omp;
  MPI_Comm comm_master;

  SubComm(comm_world, comm_omp, comm_master);

  PrintStats(comm_world, comm_omp, comm_master);

  MPI_Finalize();
}

