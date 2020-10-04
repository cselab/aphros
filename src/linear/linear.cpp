// Created by Petr Karnakov on 03.10.2020
// Copyright 2020 ETH Zurich

#include "linear.ipp"


namespace linear {

using M = MeshStructured<double, 3>;
template class SolverHypre<M>;
template class SolverConjugate<M>;
template class SolverJacobi<M>;

} // namespace linear
