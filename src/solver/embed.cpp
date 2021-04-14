// Created by Petr Karnakov on 11.02.2020
// Copyright 2020 ETH Zurich

#include "embed.ipp"

#define X(dim) template class Embed<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X
