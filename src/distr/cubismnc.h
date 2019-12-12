#pragma once

#include <mpi.h>
#include <memory>

#include "parse/vars.h"
#include "distr.h"

template <class KF>
std::unique_ptr<DistrMesh<KF>> CreateCubismnc(
    MPI_Comm comm, KF& kf, Vars& par);
