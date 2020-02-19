// Created by Petr Karnakov on 18.01.2020
// Copyright 2020 ETH Zurich

#include "proj_eb.ipp"

using M = MeshStructured<double, 3>;
using EB = Embed<M>;
template class ProjEmbed<EB>;
