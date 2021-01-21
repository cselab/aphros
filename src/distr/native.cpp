// Created by Petr Karnakov on 02.01.2021
// Copyright 2021 ETH Zurich

#include <mpi.h>
#include <sstream>

#include "geom/mesh.h"
#include "kernel/kernelmesh.h"

#include "native.ipp"

using M = MeshStructured<double, 3>;

template std::unique_ptr<DistrMesh<M>> CreateNative<M>(
    MPI_Comm, const KernelMeshFactory<M>&, Vars&);
