// Created by Petr Karnakov on 01.09.2019
// Copyright 2019 ETH Zurich

#include "vofm.ipp"
#include "embed.h"

#define X(dim) template class Vofm<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X

#define X(dim) template class Vofm<Embed<MeshCartesian<double, dim>>>;
MULTIDIMX
#undef X
