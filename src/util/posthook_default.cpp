// Created by Petr Karnakov on 25.02.2020
// Copyright 2020 ETH Zurich

#include "posthook.h"

template <class M>
void PostHook(const Vars&, const FieldCell<typename M::Vect>&, M&) {}

template <class M>
void PostHook(
    const Vars&, const FieldCell<typename M::Vect>&, M&, const Embed<M>&) {}

template <class M>
void InitVelHook(FieldCell<typename M::Vect>&, const Vars&, const M&) {}

template <class M>
void InitVelHook(
    FieldCell<typename M::Vect>&, const Vars&, const M&, const Embed<M>&) {}

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using EB = Embed<M>;

template void PostHook(const Vars&, const FieldCell<Vect>&, M&);
template void PostHook(const Vars&, const FieldCell<Vect>&, M&, const EB&);
template void InitVelHook(FieldCell<Vect>&, const Vars&, const M&);
template void InitVelHook(FieldCell<Vect>&, const Vars&, const M&, const EB&);
