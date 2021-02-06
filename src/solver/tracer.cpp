// Created by Petr Karnakov on 24.06.2020
// Copyright 2020 ETH Zurich

#include "tracer.ipp"
#include "embed.h"

#define X(dim) template class Tracer<MeshStructured<double, dim>>;
MULTIDIMX
#undef X

#define X(dim) template class Tracer<Embed<MeshStructured<double, dim>>>;
MULTIDIMX
#undef X
