// Created by Petr Karnakov on 30.05.2020
// Copyright 2020 ETH Zurich

#include <fstream>
#include <iostream>
#include <limits>

#include <func/init_u.h>
#include <func/init_vel.h>
#include <kernel/hydro.h>
#include <solver/multi.h>
#include <util/format.h>
#include <util/posthook.h>

template <class M>
void FluidDummyHook(
    FieldCell<typename M::Vect>& fcvel, FieldFace<typename M::Scal>& ffv,
    typename M::Scal t, typename M::Scal dt, const Vars& var, const M& m) {
  using Scal = typename M::Scal;
  const Scal rt = var.Double["reverse_t"];
  if (t >= rt && t < rt + dt) {
    for (auto f : m.Faces()) {
      ffv[f] *= -1;
    }
    for (auto c : m.AllCells()) {
      fcvel[c] *= -1;
    }
  }
}

template <class M>
void InitHook(Hydro<M>* hydro) {
  auto& var = hydro->var;
  auto mod = [&var](auto& fcu, auto& fccl, auto layers, auto& m) {
    auto buf = ReadPrimList(var.String["list_path"], m.IsRoot());
    InitOverlappingComponents(buf, fcu, fccl, layers, m);
  };
  if (auto* as = dynamic_cast<Vofm<M>*>(hydro->as_.get())) {
    as->AddModifier(mod);
  }
}

using M = MeshCartesian<double, 3>;

template void InitHook(Hydro<M>*);
template void FluidDummyHook(
    FieldCell<typename M::Vect>&, FieldFace<typename M::Scal>&,
    typename M::Scal t, typename M::Scal dt, const Vars&, const M&);
