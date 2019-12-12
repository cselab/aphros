#pragma once

#include <mpi.h>
#include <memory>

#include "distr.h"
#include "parse/vars.h"

template <class KF>
std::unique_ptr<DistrMesh<KF>> CreateCubism(MPI_Comm, KF&, Vars&);
