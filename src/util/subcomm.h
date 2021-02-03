// Created by Petr Karnakov on 28.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include "util/mpi.h"

// Returns communicator for hybrid shared-memory kernels.
// comm_world: copy of MPI_COMM_WORLD
// comm_omp: communicator with all ranks on the same node
// comm_master: communicator with one master rank from each node
// (see example in subcomm.cpp)
void SubComm(MPI_Comm& comm_world, MPI_Comm& comm_omp, MPI_Comm& comm_master);
void PrintStats(MPI_Comm comm_world, MPI_Comm comm_omp, MPI_Comm comm_master);
// Moves caller thread to given cpu
void SetAffinity(int cpu);
