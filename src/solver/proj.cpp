// Created by Petr Karnakov on 30.09.2019
// Copyright 2019 ETH Zurich

#include "proj.ipp"
#include "embed.h"

using M = MeshStructured<double, 3>;
template class Proj<M>;
template class Proj<Embed<M>>;
