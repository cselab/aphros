#undef NDEBUG
#include <iostream>
#include <vector>
#include <mpi.h>
#include <cassert>

#include "util/subcomm.h"

#define EV(x) (#x) << "=" << (x) << " "

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int size, rank;

  auto comm = MPI_COMM_WORLD;

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  MPI_Comm comm_world;
  MPI_Comm comm_omp;
  MPI_Comm comm_master;

  SubComm(comm_world, comm_omp, comm_master);

  std::cout
      << EV(size)
      << EV(rank)
      << std::endl;

  MPI_Finalize();
}

