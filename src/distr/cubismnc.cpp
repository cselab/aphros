// Created by Petr Karnakov on 25.11.2019
// Copyright 2019 ETH Zurich

#include <sstream>

#include "geom/mesh.h"
#include "kernel/kernelmesh.h"
#include "parse/vars.h"

#include "cubismnc.ipp"

DECLARE_FORCE_LINK_TARGET(distr_cubismnc);

template <size_t bx, size_t by, size_t bz, class M>
static void Try(
    MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var,
    std::unique_ptr<DistrMesh<M>>& r) {
  using Scal = typename M::Scal;
  if (r) return;
  // using Par0 = GPar<Scal, bx, by, bz, 0>; // 0 halos on either side
  // using Par1 = GPar<Scal, bx, by, bz, 1>; // 1 halos on either side
  using Par2 = GPar<Scal, bx, by, bz, 2>; // 2 halos on either side
  // using Par3 = GPar<Scal, bx, by, bz, 3>; // 3 halos on either side

  if (var.Int["bsx"] == bx && var.Int["bsy"] == by && var.Int["bsz"] == bz) {
    switch (var.Int["hl"]) {
      // case 0:
      //  r.reset(new Cubismnc<Par0, M>(comm, kf, var));
      // case 1:
      //  r.reset(new Cubismnc<Par1, M>(comm, kf, var));
      case 2:
        r.reset(new Cubismnc<Par2, M>(comm, kf, var));
      // case 3:
      //  r.reset(new Cubismnc<Par3, M>(comm, kf, var));
      default:
        break;
    }
  }
}

template <class M>
class ModuleDistrCubismnc : public ModuleDistr<M> {
 public:
  ModuleDistrCubismnc() : ModuleDistr<M>("cubismnc") {}
  std::unique_ptr<DistrMesh<M>> Make(
      MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var) override {
    std::unique_ptr<DistrMesh<M>> r;
    // 3D
    Try<32, 32, 32>(comm, kf, var, r);
    Try<16, 16, 16>(comm, kf, var, r);
    Try<8, 8, 8>(comm, kf, var, r);
    // 2D
    Try<32, 32, 1>(comm, kf, var, r);
    Try<16, 16, 1>(comm, kf, var, r);
    Try<8, 8, 1>(comm, kf, var, r);
    fassert(
        r, util::Format(
               "{}: no instance with bs={} {} {}, hl={}", __func__,
               var.Int["bsx"], var.Int["bsy"], var.Int["bsz"], var.Int["hl"]));
    return r;
  }
};

using M = MeshCartesian<double, 3>;

bool kRegDistrCubismnc[] = {
    RegisterModule<ModuleDistrCubismnc<M>>(),
};
