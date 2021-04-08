// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#include "convdiffvg.ipp"
#include "convdiffe.h"
#include "convdiffi.h"

#define X(dim)                        \
  template class ConvDiffVectGeneric< \
      MeshCartesian<double, dim>,    \
      ConvDiffScalImp<MeshCartesian<double, dim>>>;
MULTIDIMX
#undef X

#define X(dim)                        \
  template class ConvDiffVectGeneric< \
      MeshCartesian<double, dim>,    \
      ConvDiffScalExp<MeshCartesian<double, dim>>>;
MULTIDIMX
#undef X

#define X(dim)                            \
  template class ConvDiffVectGeneric<     \
      Embed<MeshCartesian<double, dim>>, \
      ConvDiffScalImp<Embed<MeshCartesian<double, dim>>>>;
MULTIDIMX
#undef X

#define X(dim)                            \
  template class ConvDiffVectGeneric<     \
      Embed<MeshCartesian<double, dim>>, \
      ConvDiffScalExp<Embed<MeshCartesian<double, dim>>>>;
MULTIDIMX
#undef X
