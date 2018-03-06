#pragma once

#include <array>
#include <memory>
#include <mpi.h>

#include "Kernel.h"
#include "Vars.h"
#include "Distr.h"

Distr* CreateLocal(MPI_Comm comm, KernelFactory& kf, Vars& par);
