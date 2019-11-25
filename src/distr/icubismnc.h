#pragma once

#include <mpi.h>

#include "kernel/kernel.h"
#include "parse/vars.h"
#include "distr.h"

Distr* CreateCubismnc(MPI_Comm comm, KernelFactory& kf, Vars& par);
