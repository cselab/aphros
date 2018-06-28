#pragma once

#include <limits>

#include "solver/solver.h"

namespace solver {

// Assign colors to connected sets with u > 0.
template <class M_>
class Tracker {
  using M = M_;
  using Scal = typename M::Scal;

 public:
  // fco: initial color, [g]roup [a]
  // th: threshold for detection u > 0
  // dim: space dimension, 2 or 3
  Tracker(M& m, const FieldCell<Scal>& fco, Scal th, size_t dim)
      : m(m), fco_(fco), th_(th), dim_(dim) {}
  // Propagates color.
  // fcu: volume fraction [a]
  void Update(const FieldCell<Scal>& fcu);
  // Returns color [a]
  const FieldCell<Scal>& GetColor() const { return fco_; }

 private:
  M& m;
  FieldCell<Scal> fco_; // colo[r] [a]
  size_t dim_;
  const Scal kNone = std::numeric_limits<Scal>::max(); // no color
};

template <class M_>
void Tracker<M_>::Update(const FieldCell<Scal>& fcu) {
  using MIdx = typename M::MIdx;
  auto& bc = m.GetBlockCells();

  const size_t sw = 1; // stencil halfwidth, [sw,sw]
  const int sn = sw * 2 + 1; // stencil size
  // block of offsets
  GBlock<IdxCell, dim> bo(MIdx(-sw, -sw, dim_ == 2 ? 0 : -sw), 
                          MIdx(sn, sn, dim_ == 2 ? 1 : sn)); 

  for (auto c : m.Cells()) {
    Scal& o = fco_[c];
    MIdx w = bc.GetMIdx(c);
    if (fcu[c] > th_) { // cell contains liquid
      // traverse neighbours, pick smallest color
      for (MIdx wo : bo) {
        IdxCell cn = bc.GetIdx(w + wo); // neighbour cell
        if (fcu[cn] > th) {
          o = std::min(o, fco_[cn]);
        }
      }
    } else {  // no liquid in cell
      o = kNone;
    }
  }
}


} // namespace solver
