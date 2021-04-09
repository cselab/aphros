// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include "kernelmesh.h"

template <class M_, class Par_>
class KernelMeshPar : public KernelMesh<M_> {
  using P = KernelMesh<M_>;

 public:
  using M = M_;
  using Par = Par_;

  KernelMeshPar(Vars& var_, const generic::BlockInfoProxy<M::dim>& bi, Par& par)
      : KernelMesh<M>(var_, bi), par_(par) {}
  void Run() override = 0;

 protected:
  using P::IsLead;
  using P::IsRoot;
  using P::m;
  using P::var;
  using P::var_mutable;
  Par& par_;
};

template <class M_, class Kernel_>
class KernelMeshParFactory : public KernelMeshFactory<M_> {
 public:
  using M = M_;
  using Kernel = Kernel_;
  using Par = typename Kernel::Par;
  KernelMeshParFactory(Par& par) : par_(par) {}
  Kernel* Make(
      Vars& var, const generic::BlockInfoProxy<M::dim>& p) const override {
    return new Kernel(var, p, par_);
  }

 protected:
  Par& par_;
};
