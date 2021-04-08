// Created by Petr Karnakov on 27.10.2020
// Copyright 2020 ETH Zurich

#include "mesh.ipp"

#define XX(M)                                                      \
  template class M;                                                \
  template void M::ApplyNanFaces(FieldCell<typename M::Scal>& fc); \
  template void M::ApplyNanFaces(FieldCell<typename M::Vect>& fc);

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX
