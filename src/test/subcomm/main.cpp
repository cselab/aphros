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

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  std::cout
      << EV(size)
      << EV(rank)
      << std::endl;

  MPI_Finalize();
}

