// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#include "convdiffvg.ipp"
#include "convdiffe.h"
#include "convdiffi.h"

using M = MeshStructured<double, 3>;
using EB = Embed<M>;

template class ConvDiffVectGeneric<M, ConvDiffScalImp<M>>;
template class ConvDiffVectGeneric<M, ConvDiffScalExp<M>>;
//template class ConvDiffVectGeneric<EB, ConvDiffScalImp<EB>>;
template class ConvDiffVectGeneric<EB, ConvDiffScalExp<EB>>;
