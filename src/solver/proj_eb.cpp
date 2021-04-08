// Created by Petr Karnakov on 30.09.2019
// Copyright 2019 ETH Zurich

#include "embed.h"
#include "proj.ipp"

#define X(dim) template class Proj<Embed<MeshCartesian<double, dim>>>;
MULTIDIMX
#undef X
