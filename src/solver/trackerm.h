#pragma once

#include <limits>

#include "solver/solver.h"
#include "solver/multi.h"

namespace solver {

// Assign colors to connected sets with u > 0.
template <class M_>
class Trackerm {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  static const size_t dim = M::dim;

  Trackerm(M& m, const GRange<size_t>& layers)
      : m(m), layers(layers)
      , fcclm_(m, kClNone), fcimm_(m, 0), fcim_(m, 0) {}
  void Update(const FieldCell<Scal>& fccl);
  // Returns image vector
  const FieldCell<Vect>& GetImage() const { return fcim_; }
  static constexpr Scal kNone = -1.; // no color

  struct Bit {
    int w0 : 16;
    int w1 : 16;
    int w2 : 16;
  };

  union Union {
    Scal a;
    Bit b;
  };

  static Scal Pack(MIdx w) {
    Union u;
    u.b.w0 = w[0];
    u.b.w1 = w[1];
    u.b.w2 = w[2];
    return u.a;
  }

  static MIdx Unpack(Scal a) {
    Union u;
    u.a = a;
    return MIdx(u.b.w0, u.b.w1, u.b.w2);
  }

 private:
  M& m;
  static constexpr Scal kClNone = -1;
  const GRange<size_t>& layers;
  Multi<FieldCell<Scal>> fcclm_; // color previous [a]
  Multi<FieldCell<Vect>> fcimm_; // image previous [a]
  Multi<FieldCell<Vect>> fcim_;  // image current [a]

};

template <class M_>
void Trackerm<M_>::Update(const FieldCell<Scal>& fcu) {
  //auto& bc = m.GetIndexCells();

  auto sem = m.GetSem("upd");

  if (sem("")) {
    /*
    const int sw = 1; // stencil halfwidth, [-sw,sw]
    const int sn = sw * 2 + 1; // stencil size
    GBlock<IdxCell, dim> bo(MIdx(-sw), MIdx(sn));

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
        fcim_[c] = Vect(0);
      }
    }
    m.Comm(&fcim_);
    */
  }
}

} // namespace solver
