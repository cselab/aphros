// Created by Petr Karnakov on 04.10.2020
// Copyright 2020 ETH Zurich

#include "linear.ipp"

#define X(dim) template class ULinear<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X
