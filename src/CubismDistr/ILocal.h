#pragma once

#include <array>
#include <memory>
#include <mpi.h>

#include "Kernel.h"
#include "Vars.h"

std::unique_ptr<Distr> CreateLocal(
    MPI_Comm comm, KernelFactory& kf, int bs, int es, int h, Vars& par);
