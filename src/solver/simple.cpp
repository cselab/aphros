// Created by Petr Karnakov on 31.07.2018
// Copyright 2018 ETH Zurich

#include "simple.ipp"

using M = MeshStructured<double, 3>;
template class Simple<M>;
template class Simple<Embed<M>>;
