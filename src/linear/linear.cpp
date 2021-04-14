// Created by Petr Karnakov on 03.10.2020
// Copyright 2020 ETH Zurich

#include "linear.ipp"

namespace linear {

#define X(dim) template class SolverConjugate<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X

#define X(dim) template class SolverJacobi<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X

#define X(dim) \
  RegisterModule<ModuleLinearConjugate<MeshCartesian<double, dim>>>(),
bool kReg_conjugate[] = {MULTIDIMX};
#undef X

#define X(dim) \
  RegisterModule<ModuleLinearJacobi<MeshCartesian<double, dim>>>(),
bool kReg_jacobi[] = {MULTIDIMX};
#undef X

} // namespace linear
