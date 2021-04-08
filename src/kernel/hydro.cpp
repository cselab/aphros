// Created by Petr Karnakov on 31.07.2018
// Copyright 2018 ETH Zurich

#include "hydro.ipp"

#define X(dim) template class Hydro<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X
