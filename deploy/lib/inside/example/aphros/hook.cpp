// Created by Petr Karnakov on 22.03.2020
// Copyright 2020 ETH Zurich

#include <iostream>

#include <util/posthook.h>

template <class M>
void InitEmbedHook(
    FieldNode<typename M::Scal>& fn_levelset, const Vars&, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  if (m.IsRoot()) {
    std::cout << "Embedded boundaries from lib/inside"
              << std::endl;
  }
  for (auto n : m.AllNodes()) {
    (void) n;
  }
}

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

template void InitEmbedHook(FieldNode<Scal>&, const Vars&, const M&);
