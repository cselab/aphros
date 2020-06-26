// Created by Petr Karnakov on 26.06.2020
// Copyright 2020 ETH Zurich

#include "particles.ipp"
#include "embed.h"

using M = MeshStructured<double, 3>;
template class Particles<M>;
template class Particles<Embed<M>>;
