// Created by Petr Karnakov on 20.03.2021
// Copyright 2021 ETH Zurich

#include "hydro_post.ipp"

#define X(dim) template struct HydroPost<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X
