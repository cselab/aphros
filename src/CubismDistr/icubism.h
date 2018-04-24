#pragma once

#include <array>
#include <mpi.h>
#include <memory>

#include "kernel.h"
#include "vars.h"
#include "distr.h"

Distr* CreateCubism(MPI_Comm comm, KernelFactory& kf, Vars& par);
