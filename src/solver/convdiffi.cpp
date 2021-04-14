// Created by Petr Karnakov on 30.07.2018
// Copyright 2018 ETH Zurich

#include "convdiffi.ipp"

#define X(dim) template class ConvDiffScalImp<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X

#define X(dim) \
  template class ConvDiffScalImp<Embed<MeshCartesian<double, dim>>>;
MULTIDIMX
#undef X
