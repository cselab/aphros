// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#include "convdiffvg.ipp"
#include "convdiffe.h"
#include "convdiffi.h"

template class ConvDiffVectGeneric<
    MeshStructured<double, 3>, ConvDiffScalImp<MeshStructured<double, 3>>>;

template class ConvDiffVectGeneric<
    MeshStructured<double, 3>, ConvDiffScalExp<MeshStructured<double, 3>>>;
