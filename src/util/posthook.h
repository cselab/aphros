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
