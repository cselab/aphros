#pragma once

#include <memory>

#include <geom/mesh.h>

namespace solver {

template <class M_>
class UNormal {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  // Computes normal by combined Young's scheme and height-functions
  // m: mesh
  // fcu: volume fraction
  // fci: interface mask (1: contains interface)
  // edim: effective dimension
  // Output: set to NaN if fci=0
  // fcn: normal with norm1()=1, antigradient of fcu [s] 
  // fck: curvature
  static void CalcNormal(
      M& m, const FieldCell<Scal>& fcu, const FieldCell<bool>& fci,
      size_t edim, FieldCell<Vect>& fcn, FieldCell<Scal>& fck);

 private:
  struct Imp; // implementation
};

} // namespace solver
