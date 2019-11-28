#undef NDEBUG
#include <iostream>
#include <vector>
#include <mpi.h>
#include <cassert>

#include "util/subcomm.h"

#define EV(x) (#x) << "=" << (x) << " "

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int size_world, rank_world;
  int size_omp, rank_omp;
  int size_master, rank_master;

  MPI_Comm comm_world;
  MPI_Comm comm_omp;
  MPI_Comm comm_master;

  SubComm(comm_world, comm_omp, comm_master);

  MPI_Comm_size(comm_world, &size_world);
  MPI_Comm_rank(comm_world, &rank_world);
  MPI_Comm_size(comm_omp, &size_omp);
  MPI_Comm_rank(comm_omp, &rank_omp);

  if (rank_omp == 0) {
    MPI_Comm_size(comm_master, &size_master);
    MPI_Comm_rank(comm_master, &rank_master);
  }

  std::cout
      << EV(size_world) << EV(rank_world) << std::endl
      << EV(size_omp) << EV(rank_omp) << std::endl
      << EV(size_master) << EV(rank_master) << std::endl
      << std::endl;

  MPI_Finalize();
}

