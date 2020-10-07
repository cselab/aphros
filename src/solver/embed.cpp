// Created by Petr Karnakov on 11.02.2020
// Copyright 2020 ETH Zurich

#include "embed.ipp"

using M = MeshStructured<double, 3>;
using EB = Embed<M>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

template class Embed<M>;
