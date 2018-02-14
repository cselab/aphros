#pragma once

#include <array>
#include <mpi.h>
#include <memory>

#include "Kernel.h"

std::unique_ptr<Distr> CreateCubism(
    MPI_Comm comm, KernelFactory& kf, 
    int bs, Idx b, Idx p, int es, int h);
