#pragma once

#include <limits>

#include "solver/solver.h"

namespace solver {

// Embedded boundaries.
template <class M_>
class Embed {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;

 public:
  // fct: cell type (0: regular, 1: cut, 2: excluded)
  // fcn: normal
  // fca: plane constant
  Embed(M& m, const FieldCell<char>& fct, 
        const FieldCell<Vect>& fcn, const FieldCell<Scal>& fca)
      : m(m), fct_(fct), fcn_(fcn), fca_(fca) {}
  const FieldCell<Scal>& GetCellType() const { return fct_; }
  const FieldCell<Scal>& GetNormal() const { return fcn_; }
  const FieldCell<Vect>& GetPlane() const { return fca_; }
  static constexpr Scal kNone = -1.; // no color

 private:
  M& m;
  FieldCell<char> fct_; // cell type [a] (0: regular, 1: cut, 2: excluded)
  FieldCell<Vect> fcn_; // normal [a]
  FieldCell<Scal> fca_; // plane constant [a]
};


} // namespace solver
