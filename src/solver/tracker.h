#pragma once

#include <limits>

#include "solver/solver.h"

namespace solver {

// Assign colors to connected sets with u > 0.
template <class M_>
class Tracker {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;

 public:
  // fccl: initial color, [g]roup [a]
  // th: threshold for detection u > 0
  // dm: space dimension, 2 or 3
  Tracker(M& m, const FieldCell<Scal>& fccl, Scal th, size_t dm)
      : m(m), fccl_(fccl), fcim_(m, Vect(0)), th_(th), dm_(dm) {}
  // Propagates color.
  // fcu: volume fraction [a]
  void Update(const FieldCell<Scal>& fcu);
  // Returns color [a]
  const FieldCell<Scal>& GetColor() const { return fccl_; }
  const FieldCell<Vect>& GetImage() const { return fcim_; }
  static constexpr Scal kNone = -1.; // no color

 private:
  M& m;
  FieldCell<Scal> fccl_; // color [a]
  FieldCell<Vect> fcim_; // image  [a]
  size_t dm_;
  Scal th_;
};

template <class M_>
void Tracker<M_>::Update(const FieldCell<Scal>& fcu) {
  using MIdx = typename M::MIdx;
  auto& bc = m.GetBlockCells();

  auto sem = m.GetSem("upd");

  if (sem("")) {
    const size_t sw = 1; // stencil halfwidth, [sw,sw]
    const int sn = sw * 2 + 1; // stencil size
    // block of offsets
    GBlock<IdxCell, dim> bo(MIdx(-sw, -sw, dm_ == 2 ? 0 : -sw), 
                            MIdx(sn, sn, dm_ == 2 ? 1 : sn)); 

    MIdx gs = m.GetGlobalSize();
    for (auto c : m.Cells()) {
      Scal& o = fccl_[c];
      MIdx w = bc.GetMIdx(c);
      if (fcu[c] > th_) { // liquid in cell
        // traverse neighbours, pick minimal color
        for (MIdx wo : bo) {
          MIdx wn = w + wo;
          IdxCell cn = bc.GetIdx(w + wo); // neighbour 
          Scal on = fccl_[cn];
          if (fcu[cn] > th_ && on != kNone) { // liquid and color in neighbour
            if (o == kNone || on < o) { // update if empty or smaller
              o = on;
              Vect im = fcim_[cn];
              for (size_t d = 0; d < dim; ++d) { // periodic
                (wn[d] < 0) && (im[d] += 1.);
                (wn[d] >= gs[d]) && (im[d] -= 1.);
              }
              fcim_[c] = im; // image from neighbour
            }
          }
        }
      } else {  // no liquid in cell
        o = kNone;
      }
    }
    m.Comm(&fccl_);
  }
}


} // namespace solver
