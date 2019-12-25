#pragma once

#include <memory>

#include "geom/mesh.h"

template <class M_>
struct UCurv {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  // Computes curvature with height functions.
  // fcu: volume fraction [a]
  // fcn: normal
  // edim: effective dimension (2 or 3)
  // Output:
  // fck: curvature
  static void CalcCurvHeight(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcn, size_t edim,
      FieldCell<Scal>& fck, M& m);

 private:
  struct Imp;
};
