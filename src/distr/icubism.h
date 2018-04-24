#pragma once

#include <array>
#include <mpi.h>
#include <memory>

#include "kernel/kernel.h"
#include "parse/vars.h"
#include "distr.h"

Distr* CreateCubism(MPI_Comm comm, KernelFactory& kf, Vars& par);
