// Created by Petr Karnakov on 25.02.2020
// Copyright 2020 ETH Zurich

#include <iostream>

#include <util/posthook.h>

template <class M>
void PostHook(
    const Vars&, const FieldCell<typename M::Vect>&, M&, const Embed<M>&) {
  std::cout << "USER_DEFINED_POSTHOOK" << std::endl;
}

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using EB = Embed<M>;

template void PostHook(const Vars&, const FieldCell<Vect>&, M&, const EB&);
