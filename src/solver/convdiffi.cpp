// Created by Petr Karnakov on 30.07.2018
// Copyright 2018 ETH Zurich

#include "convdiffi.ipp"

using M = MeshStructured<double, 3>;
template class ConvDiffScalImp<M>;
//template class ConvDiffScalImp<Embed<M>>;
