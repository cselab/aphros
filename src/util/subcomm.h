#pragma once

#include <mpi.h>

// Returns communicator for OpenMP kernels.
// comm_world: copy of MPI_COMM_WORLD
// comm_omp: communicator with all ranks on the same node
// comm_master: communicator wiht one master rank from each node
// (see example in subcomm.cpp)
void SplitComm(
    MPI_Comm& comm_world, MPI_Comm& comm_omp, MPI_Comm& comm_master);
