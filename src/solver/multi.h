#pragma once

#include <limits>

#include "solver/solver.h"

namespace solver {

// Assign colors to connected sets with u > 0.
template <class M_>
class Multi {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;

 public:
  // fccl: color [a]
  Multi(M& m, const FieldCell<Scal>* fccl, size_t dm) 
      : m(m), fccl_(fccl), fcmask_(m, 0), dm_(dm) {}
  // Propagates color.
  // fcu: volume fraction [a]
  void Update(const FieldCell<Scal>& fcu);
  // Returns color [a]
  const FieldCell<Scal>& GetMask() const { return fcmask_; }
  static constexpr Scal kNone = -1.; // no color

 private:
  M& m;
  const FieldCell<Scal>* fccl_; // color
  FieldCell<Scal> fcmask_; // mask
  size_t dm_;
};

template <class M_>
void Multi<M_>::Update(const FieldCell<Scal>& fcu) {
  using MIdx = typename M::MIdx;
  auto& bc = m.GetIndexCells();

  auto sem = m.GetSem("upd");

  if (sem("")) {
    const int sw = 1; // stencil halfwidth, [sw,sw]
    const int sn = sw * 2 + 1; // stencil size
    // block of offsets
    GBlock<IdxCell, dim> bo(MIdx(-sw, -sw, dm_ == 2 ? 0 : -sw), 
                            MIdx(sn, sn, dm_ == 2 ? 1 : sn)); 
    (void) bo;
    fcmask_ = fcu;
  }
}


} // namespace solver
