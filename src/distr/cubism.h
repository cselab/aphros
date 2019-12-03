#pragma once

#include <mpi.h>
#include <memory>

#include "kernel/kernel.h"
#include "parse/vars.h"
#include "distr.h"

std::unique_ptr<Distr> CreateCubism(
    MPI_Comm comm, KernelFactory& kf, Vars& par);
