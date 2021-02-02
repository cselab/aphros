// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#include <mpi.h>
#include <sstream>

#include "geom/mesh.h"
#include "kernel/kernelmesh.h"
#include "parse/vars.h"

#include "local.ipp"

DECLARE_FORCE_LINK_TARGET(distr_local);

template <class M>
class ModuleDistrLocal : public ModuleDistr<M> {
 public:
  ModuleDistrLocal() : ModuleDistr<M>("local") {}
  std::unique_ptr<DistrMesh<M>> Make(
      MPI_Comm comm, const KernelMeshFactory<M>& factory, Vars& var) override {
    return std::unique_ptr<DistrMesh<M>>(new Local<M>(comm, factory, var));
  }
};

using M = MeshStructured<double, 3>;

bool kRegDistrLocal[] = {
    RegisterModule<ModuleDistrLocal<M>>(),
};
