// Created by Petr Karnakov on 30.09.2019
// Copyright 2019 ETH Zurich

#include "embed.h"
#include "proj.ipp"

using M = MeshStructured<double, 3>;
template class Proj<Embed<M>>;
