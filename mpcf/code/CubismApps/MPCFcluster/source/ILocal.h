#pragma once

#include <array>
#include <mpi.h>
#include <memory>

#include "Kernel.h"

using Idx = std::array<int, 3>;

std::unique_ptr<Distr> CreateLocal(
    MPI_Comm comm, KernelFactory& kf, 
    int bs, Idx b, Idx p, int es, int h);
