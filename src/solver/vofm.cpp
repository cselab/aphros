// Created by Petr Karnakov on 01.09.2019
// Copyright 2019 ETH Zurich

#include "vofm.ipp"
#include "embed.h"

template class Vofm<MeshStructured<double, 3>>;
template class Vofm<Embed<MeshStructured<double, 3>>>;
