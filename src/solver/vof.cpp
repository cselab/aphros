// Created by Petr Karnakov on 30.07.2018
// Copyright 2018 ETH Zurich

#include "vof.ipp"
#include "embed.h"

using M = MeshStructured<double, 3>;
template class Vof<M>;
template class Vof<Embed<M>>;
