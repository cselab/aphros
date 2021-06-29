// Created by Petr Karnakov on 25.12.2019
// Copyright 2019 ETH Zurich

#include "curv.ipp"

#define X(dim) template struct UCurv<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X

namespace curvature {

#define X(dim) template class Particles<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X

#define XX(M)                                              \
  template std::unique_ptr<Estimator<M>> MakeEstimator<M>( \
      const Vars&, M&, const GRange<size_t>&);

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX
#undef X
#undef XX

} // namespace curvature
