#pragma once

#include <mpi.h>
#include <memory>

#include "distr.h"
#include "parse/vars.h"

template <class M>
std::unique_ptr<DistrMesh<M>> CreateCubism(
    MPI_Comm, const KernelMeshFactory<M>&, Vars&);
