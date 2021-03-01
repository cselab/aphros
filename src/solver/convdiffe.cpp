// Created by Petr Karnakov on 30.04.2019
// Copyright 2019 ETH Zurich

#include "convdiffe.ipp"

#define X(dim) template class ConvDiffScalExp<MeshStructured<double, dim>>;
MULTIDIMX
#undef X

#define X(dim) \
  template class ConvDiffScalExp<Embed<MeshStructured<double, dim>>>;
MULTIDIMX
#undef X
