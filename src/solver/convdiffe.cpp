// Created by Petr Karnakov on 30.04.2019
// Copyright 2019 ETH Zurich

#include "convdiffe.ipp"

using M = MeshStructured<double, 3>;
template class ConvDiffScalExp<M>;
template class ConvDiffScalExp<Embed<M>>;
