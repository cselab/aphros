// Created by Petr Karnakov on 25.12.2019
// Copyright 2019 ETH Zurich

#include "curv.ipp"

#define X(dim) template class UCurv<MeshStructured<double, dim>>;
MULTIDIMX
#undef X
