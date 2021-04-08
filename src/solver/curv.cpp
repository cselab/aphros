// Created by Petr Karnakov on 25.12.2019
// Copyright 2019 ETH Zurich

#include "curv.ipp"

#define X(dim) template struct UCurv<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X
