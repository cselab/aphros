// Created by Petr Karnakov on 11.02.2020
// Copyright 2020 ETH Zurich

#include "approx_eb.ipp"
#include "linear/linear.h"

using M = MeshStructured<double, 3>;
using EB = Embed<M>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using T = Scal;
using TV = Vect;
constexpr size_t dim = M::dim;

template struct ULinear<Scal>;
template struct UEmbed<M>;
