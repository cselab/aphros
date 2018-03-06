#pragma once

#include <array>
#include <mpi.h>
#include <memory>

#include "Kernel.h"
#include "Vars.h"
#include "Distr.h"

Distr* CreateCubism(
    MPI_Comm comm, KernelFactory& kf, int bs, int es, int h, Vars& par);
