// Created by Petr Karnakov on 27.09.2020
// Copyright 2020 ETH Zurich

#include "electro.ipp"
#include "embed.h"

#define X(dim) template class Electro<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X

#define X(dim) template class Electro<Embed<MeshCartesian<double, dim>>>;
MULTIDIMX
#undef X
