// Created by Petr Karnakov on 02.09.2019
// Copyright 2019 ETH Zurich

#include "vof.ipp"

#define X(dim) template class UVof<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X

#define X(dim) template class ModuleLabeling<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X

#define X(dim) \
  RegisterModule<ModuleLabelingPropagation<MeshCartesian<double, dim>>>(),
bool kReg_propagation[] = {MULTIDIMX};
#undef X

#define X(dim) \
  RegisterModule<ModuleLabelingUnionFind<MeshCartesian<double, dim>>>(),
bool kReg_unionfind[] = {MULTIDIMX};
#undef X
