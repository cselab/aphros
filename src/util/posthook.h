// Created by Petr Karnakov on 25.02.2020
// Copyright 2020 ETH Zurich

#pragma once

#include "geom/mesh.h"
#include "parse/vars.h"
#include "solver/embed.h"

template <class M>
void PostHook(const Vars& var, const FieldCell<typename M::Vect>& fcvel, M& m);

template <class M>
void PostHook(
    const Vars& var, const FieldCell<typename M::Vect>& fcvel, M& m,
    const Embed<M>& eb);

template <class M>
void InitVelHook(
    FieldCell<typename M::Vect>& fcvel, const Vars& var, const M& m);

template <class M>
void InitVelHook(
    FieldCell<typename M::Vect>& fcvel, const Vars& var, const M& m,
    const Embed<M>& eb);

template <class M>
void InitEmbedHook(
    FieldNode<typename M::Scal>& fn_levelset, const Vars& var, const M& m);

template <class M>
void FluidDummyHook(
    FieldCell<typename M::Vect>& fcvel, FieldFace<typename M::Scal>& ffv,
    typename M::Scal t, typename M::Scal dt, const Vars& var, const M& m);

template <class M>
class Hydro;

template <class M>
void StepHook(Hydro<M>*);
