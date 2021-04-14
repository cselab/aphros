// Created by Petr Karnakov on 31.05.2020
// Copyright 2020 ETH Zurich

#include "fluid_dummy.ipp"
#include "embed.h"

#define X(dim) template class FluidDummy<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X

#define X(dim) template class FluidDummy<Embed<MeshCartesian<double, dim>>>;
MULTIDIMX
#undef X
