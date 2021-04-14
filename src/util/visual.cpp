// Created by Petr Karnakov on 19.03.2021
// Copyright 2021 ETH Zurich

#include "visual.ipp"

namespace util {

#define X(dim) template struct Visual<MeshCartesian<double, dim>>;
MULTIDIMX

} // namespace util
