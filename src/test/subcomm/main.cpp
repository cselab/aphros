// Created by Petr Karnakov on 28.11.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <iostream>
#include <vector>

#include "util/mpi.h"
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
