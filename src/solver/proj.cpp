// Created by Petr Karnakov on 30.09.2019
// Copyright 2019 ETH Zurich

#include "proj.ipp"
#include "embed.h"

#define X(dim) template class Proj<MeshStructured<double, dim>>;
MULTIDIMX
#undef X
