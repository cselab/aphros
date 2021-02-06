// Created by Petr Karnakov on 28.08.2019
// Copyright 2019 ETH Zurich

#include "partstrmeshm.ipp"
#include "embed.h"

#define XX(M)                                                                \
  template class PartStrMeshM<M>;                                            \
  template void PartStrMeshM<M>::Part(const Plic& plic, const Embed<M>& eb); \
  template void PartStrMeshM<M>::Part(const Plic& plic, const M& eb);

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX
