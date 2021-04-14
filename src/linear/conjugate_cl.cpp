// Created by Petr Karnakov on 14.04.2021
// Copyright 2021 ETH Zurich

#include "conjugate_cl.ipp"

namespace linear {

#define X(dim) template class SolverConjugateCL<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X

#define X(dim) \
  RegisterModule<ModuleLinearConjugateCL<MeshCartesian<double, dim>>>(),
bool kReg_conjugate_cl[] = {MULTIDIMX};
#undef X

} // namespace linear
