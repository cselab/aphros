// Created by Petr Karnakov on 28.08.2021
// Copyright 2021 ETH Zurich

#include "util/macros.h"

#include "distr_particles.ipp"

#define XX(M)                                                       \
  template void TransferParticles(                                  \
      const std::vector<std::unique_ptr<KernelMesh<M>>>&, MPI_Comm, \
      const std::function<int(int)>&);                              \
  \

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX
#undef XX
#undef X
