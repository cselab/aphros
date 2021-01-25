// Created by Petr Karnakov on 27.09.2020
// Copyright 2020 ETH Zurich

#include "electro.ipp"
#include "embed.h"

using M = MeshStructured<double, 3>;
template class Electro<M>;
template class Electro<Embed<M>>;
