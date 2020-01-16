// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#include <mpi.h>
#include <sstream>

#include "geom/mesh.h"
#include "kernel/kernelmesh.h"
#include "parse/vars.h"

#include "local.ipp"

template <class M>
std::unique_ptr<DistrMesh<M>> CreateLocal(
    MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var) {
  return std::unique_ptr<DistrMesh<M>>(new Local<M>(comm, kf, var));
}

using M = MeshStructured<double, 3>;

template std::unique_ptr<DistrMesh<M>> CreateLocal<M>(
    MPI_Comm, const KernelMeshFactory<M>&, Vars&);
