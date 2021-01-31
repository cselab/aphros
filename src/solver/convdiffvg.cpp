// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#include "convdiffvg.ipp"
#include "convdiffe.h"
#include "convdiffi.h"

#define X(dim)                        \
  template class ConvDiffVectGeneric< \
      MeshStructured<double, dim>,    \
      ConvDiffScalImp<MeshStructured<double, dim>>>;
MULTIDIMX
#undef X

#define X(dim)                        \
  template class ConvDiffVectGeneric< \
      MeshStructured<double, dim>,    \
      ConvDiffScalExp<MeshStructured<double, dim>>>;
MULTIDIMX
#undef X

#define X(dim)                            \
  template class ConvDiffVectGeneric<     \
      Embed<MeshStructured<double, dim>>, \
      ConvDiffScalImp<Embed<MeshStructured<double, dim>>>>;
MULTIDIMX
#undef X

#define X(dim)                            \
  template class ConvDiffVectGeneric<     \
      Embed<MeshStructured<double, dim>>, \
      ConvDiffScalExp<Embed<MeshStructured<double, dim>>>>;
MULTIDIMX
#undef X
