// Created by Petr Karnakov on 26.06.2020
// Copyright 2020 ETH Zurich

#include "particles.ipp"
#include "embed.h"

#define X(dim) template class Particles<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X

#define X(dim) template class Particles<Embed<MeshCartesian<double, dim>>>;
MULTIDIMX
#undef X
