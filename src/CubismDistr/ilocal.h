#pragma once

#include <array>
#include <memory>
#include <mpi.h>

#include "kernel.h"
#include "vars.h"
#include "distr.h"

Distr* CreateLocal(MPI_Comm comm, KernelFactory& kf, Vars& par);
