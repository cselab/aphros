// Created by Petr Karnakov on 17.01.2020
// Copyright 2020 ETH Zurich

#include "convdiffv_eb.ipp"

using M = MeshStructured<double, 3>;
template class ConvDiffVectEmbed<Embed<M>>;
