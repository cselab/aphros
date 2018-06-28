#pragma once

#include <limits>

#include "solver/solver.h"

namespace solver {

// Assign colors to connected sets with u > 0.
template <class M_>
class Tracker {
  using M = M_;
  using Scal = typename M::Scal;
  using M::dim;

 public:
  // fccl: initial color, [g]roup [a]
  // th: threshold for detection u > 0
  // dm: space dimension, 2 or 3
  Tracker(M& m, const FieldCell<Scal>& fccl, Scal th, size_t dm)
      : m(m), fccl_(fccl), th_(th), dm_(dm) {}
  // Propagates color.
  // fcu: volume fraction [a]
  void Update(const FieldCell<Scal>& fcu);
  // Returns color [a]
  const FieldCell<Scal>& GetColor() const { return fccl_; }

 private:
  M& m;
  FieldCell<Scal> fccl_; // color [a]
  size_t dm_;
  Scal th_;
  const Scal kNone = std::numeric_limits<Scal>::max(); // no color
};

template <class M_>
void Tracker<M_>::Update(const FieldCell<Scal>& fcu) {
  using MIdx = typename M::MIdx;
  auto& bc = m.GetBlockCells();

  const size_t sw = 1; // stencil halfwidth, [sw,sw]
  const int sn = sw * 2 + 1; // stencil size
  // block of offsets
  GBlock<IdxCell, dim> bo(MIdx(-sw, -sw, dm_ == 2 ? 0 : -sw), 
                          MIdx(sn, sn, dm_ == 2 ? 1 : sn)); 

  for (auto c : m.Cells()) {
    Scal& o = fccl_[c];
    MIdx w = bc.GetMIdx(c);
    if (fcu[c] > th_) { // cell contains liquid
      // traverse neighbours, pick smallest color
      for (MIdx wo : bo) {
        IdxCell cn = bc.GetIdx(w + wo); // neighbour cell
        if (fcu[cn] > th_) { // if neighbout contains liquid
          o = std::min(o, fccl_[cn]);
        }
      }
    } else {  // no liquid in cell
      o = kNone;
    }
  }
}


} // namespace solver
