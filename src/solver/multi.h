#pragma once

#include <limits>
#include <memory>

#include "solver/solver.h"

namespace solver {

// Assign colors to connected sets with u > 0.
template <class M_>
class MultiMask {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;

 public:
  // fccl: color [a]
  MultiMask(M& m, const FieldCell<Scal>* fccl, size_t edim);
  ~MultiMask();
  // Propagates color.
  // fcu: volume fraction [a]
  void Update(const FieldCell<Scal>& fcu);
  // Returns color [a]
  const FieldCell<Scal>& GetMask() const;
  const FieldCell<Scal>& GetMask2() const;
  static constexpr Scal kNone = -1.; // no color

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;

};

/*
template <class T>
class Multi<>
*/


} // namespace solver
