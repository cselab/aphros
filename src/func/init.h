// Created by Petr Karnakov on 28.12.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <functional>

#include "geom/field.h"
#include "parse/vars.h"

template <class M>
void InitVf(
    FieldCell<typename M::Scal>& fcu, const Vars& var, M& m, bool verbose);

template <class M>
std::function<void(
    FieldCell<typename M::Scal>&, const FieldCell<typename M::Scal>&, const M&)>
CreateInitCl(const Vars& var, bool verb);

template <class M>
std::function<void(FieldCell<typename M::Scal>&, const M&)> CreateInitU(
    const Vars& var, bool verb);

template <class M>
std::function<void(FieldCell<typename M::Scal>&, const M&)> CreateInitSig(
    const Vars& var);
