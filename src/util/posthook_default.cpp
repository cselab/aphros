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

template <class M>
void InitEmbedHook(FieldNode<typename M::Scal>&, const Vars&, const M&) {}

template <class M>
void FluidDummyHook(
    FieldCell<typename M::Vect>&, FieldFace<typename M::Scal>&,
    typename M::Scal, typename M::Scal, const Vars&, const M&) {}

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using EB = Embed<M>;

template void PostHook(const Vars&, const FieldCell<Vect>&, M&);
template void PostHook(const Vars&, const FieldCell<Vect>&, M&, const EB&);
template void InitVelHook(FieldCell<Vect>&, const Vars&, const M&);
template void InitVelHook(FieldCell<Vect>&, const Vars&, const M&, const EB&);
template void InitEmbedHook(FieldNode<Scal>&, const Vars&, const M&);
template void FluidDummyHook(
    FieldCell<Vect>&, FieldFace<Scal>&, Scal t, Scal dt, const Vars&, const M&);
