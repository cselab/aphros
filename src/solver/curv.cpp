// Created by Petr Karnakov on 25.12.2019
// Copyright 2019 ETH Zurich

#include "curv.ipp"

namespace curvature {

#define XX(M)                                              \
  template class Particles<M>;                             \
  template class Heights<M>;                               \
  template class Hybrid<M>;                                \
  template std::unique_ptr<Estimator<M>> MakeEstimator<M>( \
      const Vars&, M&, const GRange<size_t>&);

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX
#undef X
#undef XX

} // namespace curvature
