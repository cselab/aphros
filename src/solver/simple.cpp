// Created by Petr Karnakov on 31.07.2018
// Copyright 2018 ETH Zurich

#include "simple.ipp"

#define X(dim) template class Simple<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X

#define X(dim) template class Simple<Embed<MeshCartesian<double, dim>>>;
MULTIDIMX
#undef X
