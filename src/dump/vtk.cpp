// Created by Petr Karnakov on 17.06.2022
// Copyright 2022 ETH Zurich

#include "vtk.ipp"

namespace dump {

#define X(dim) template class Vtk<generic::Vect<double, dim>>;
X(1)
X(2)
X(3)
X(4)
#undef X

} // namespace dump
