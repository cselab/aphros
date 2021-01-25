// Created by Petr Karnakov on 24.06.2020
// Copyright 2020 ETH Zurich

#include "tracer.ipp"
#include "embed.h"

using M = MeshStructured<double, 3>;
template class Tracer<M>;
template class Tracer<Embed<M>>;
