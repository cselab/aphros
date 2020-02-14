// Created by Petr Karnakov on 23.01.2020
// Copyright 2020 ETH Zurich

#include "vof_eb.ipp"
#include "vof.ipp"

template class VofEmbed<MeshStructured<double, 3>>;
template class Vof<Embed<MeshStructured<double, 3>>>;
