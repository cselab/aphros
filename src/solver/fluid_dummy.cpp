// Created by Petr Karnakov on 31.05.2020
// Copyright 2020 ETH Zurich

#include "fluid_dummy.ipp"
#include "embed.h"

using M = MeshStructured<double, 3>;
template class FluidDummy<M>;
template class FluidDummy<Embed<M>>;
