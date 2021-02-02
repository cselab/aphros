// Created by Petr Karnakov on 02.01.2021
// Copyright 2021 ETH Zurich

#include <mpi.h>
#include <sstream>

#include "geom/mesh.h"
#include "kernel/kernelmesh.h"

#include "native.ipp"

DECLARE_FORCE_LINK_TARGET(distr_native);

template <class M>
class ModuleDistrNative : public ModuleDistr<M> {
 public:
  ModuleDistrNative() : ModuleDistr<M>("native") {}
  std::unique_ptr<DistrMesh<M>> Make(
      MPI_Comm comm, const KernelMeshFactory<M>& factory, Vars& var) override {
    return std::unique_ptr<DistrMesh<M>>(new Native<M>(comm, factory, var));
  }
};

using M = MeshStructured<double, 3>;

bool kRegDistrNative[] = {
    RegisterModule<ModuleDistrNative<M>>(),
};
