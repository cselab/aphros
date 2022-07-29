// Created by Petr Karnakov on 28.12.2019
// Copyright 2019 ETH Zurich

#include "init.h"
#include "init.ipp"

#define XX(M)                                                                 \
  template void InitVf(                                                       \
      FieldCell<typename M::Scal>& fcu, const Vars& var, M& m, bool verbose); \
  template std::function<void(                                                \
      FieldCell<typename M::Scal>&, const FieldCell<typename M::Scal>&,       \
      const M&)>                                                              \
  CreateInitCl(const Vars& par, bool verb);                                   \
  template std::function<void(FieldCell<typename M::Scal>&, const M&)>        \
  CreateInitU(const Vars& par, bool verb);                                    \
  template std::function<void(FieldCell<typename M::Scal>&, const M&)>        \
  CreateInitSig(const Vars& var);                                             \
  template void InitOverlappingComponents(                                    \
      std::istream& primlist, const Multi<FieldCell<typename M::Scal>*>& fcu, \
      const Multi<FieldCell<typename M::Scal>*>& fccl,                        \
      const GRange<size_t>& layers, const M& m);
#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX
