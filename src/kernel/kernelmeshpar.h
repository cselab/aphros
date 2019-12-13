#pragma once

#include "kernelmesh.h"

// Kernel aware of mesh with structure Par_ passed to constructor.
template <class M_, class Par_>
class KernelMeshPar : public KernelMesh<M_> {
 public:
  using P = KernelMesh<M_>; // parent
  using M = M_;
  using Par = Par_;
  static constexpr size_t dim = M::dim;

  KernelMeshPar(Vars& var, const MyBlockInfo& bi, Par& par)
      : KernelMesh<M>(var, bi), par_(par) {}
  void Run() override = 0;

 protected:
  using P::IsLead;
  using P::IsRoot;
  using P::m;
  using P::var;
  using P::var_mutable;
  Par& par_;
};

// Factory for KernelMeshPar.
// M_: mesh
// K_: kernel derived from KernelMeshPar with defined Par
template <class M_, class K_>
class KernelMeshParFactory : public KernelMeshFactory<M_> {
 public:
  using M = M_;
  using K = K_;
  using Par = typename K::Par;
  KernelMeshParFactory(Par& par) : par_(par) {}
  K* Make(Vars& var, const MyBlockInfo& bi) const override {
    return new K(var, bi, par_);
  }

 protected:
  Par& par_;
};
