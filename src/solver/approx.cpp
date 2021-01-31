// Created by Petr Karnakov on 29.12.2019
// Copyright 2019 ETH Zurich

#include "approx.ipp"
#include "linear/linear.h"

using Scal = double;

template std::vector<Scal> GetGradCoeffs(Scal x, const std::vector<Scal>& z);
template std::vector<Scal> GetGradCoeffs(
    Scal x, const std::vector<Scal>& z, size_t b);
template std::array<Scal, 3> GetCoeff(ConvSc);

#define XX(M)                                                   \
  template void SmoothenNode(                                   \
      FieldCell<typename M::Scal>& fc, M& m, size_t rep);       \
  template void SmoothenNode(                                   \
      FieldNode<typename M::Scal>& fn, M& m, size_t iters);     \
                                                                \
  template void BcApply(                                        \
      FieldCell<typename M::Scal>& uc,                          \
      const MapEmbed<BCond<typename M::Scal>>& me, const M& m); \
  template void BcApply(                                        \
      FieldCell<typename M::Vect>& uc,                          \
      const MapEmbed<BCond<typename M::Vect>>& me, const M& m); \
                                                                \
  template void BcReflectAll(                                   \
      FieldCell<typename M::Scal>& uc,                          \
      const MapEmbed<BCond<typename M::Scal>>& me, const M& m);

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX
