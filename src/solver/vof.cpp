// Created by Petr Karnakov on 30.07.2018
// Copyright 2018 ETH Zurich

#include "vof.ipp"
#include "embed.h"

#define X(dim) template class Vof<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X

#define X(dim) template class Vof<Embed<MeshCartesian<double, dim>>>;
MULTIDIMX
#undef X
