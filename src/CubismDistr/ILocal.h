#pragma once

#include <array>
#include <memory>
#include <mpi.h>

#include "Kernel.h"

std::unique_ptr<Distr> CreateLocal(
    MPI_Comm comm, KernelFactory& kf, int bs, Idx b, Idx p, int es, int h);
