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
      : m(m), layers(layers), fcclm_(layers), fcimm_(layers), fcim_(layers) {
    fcclm_.InitAll(FieldCell<Scal>(m, kClNone));
    fcimm_.InitAll(FieldCell<Scal>(m, Pack(MIdx(0))));
    fcim_.InitAll(FieldCell<Scal>(m, Pack(MIdx(0))));
  }
  void Update(const Multi<const FieldCell<Scal>*>& fccl);
  // Returns image vector
  Multi<const FieldCell<Scal>*> GetImage() const {
    Multi<const FieldCell<Scal>*> r(layers);
    for (auto l : layers) {
      r[l] = &fcim_[l];
    }
    return r;
  }
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
  Multi<FieldCell<Scal>> fcimm_; // image previous [a]
  Multi<FieldCell<Scal>> fcim_;  // image current [a]

};

template <class M_>
void Trackerm<M_>::Update(const Multi<const FieldCell<Scal>*>& fccl) {
  auto sem = m.GetSem("trackerm");
  if (sem("update")) {
    MIdx gs = m.GetGlobalSize();
    auto& bc = m.GetIndexCells();

    for (auto c : m.Cells()) {
      for (auto l : layers) { // check if new color appeared
        if ((*fccl[l])[c] == kClNone) continue;
        bool fndm = false;
        for (auto lm : layers) {
          if (fcclm_[lm][c] == (*fccl[l])[c]) {
            fndm = true;
            break;
          }
        }
        if (!fndm) { // new color, find same color in neighbors
          bool fndn = false;
          for (auto q : m.Nci(c)) {
            auto cn = m.GetCell(c, q);
            for (auto ln : layers) {
              if (fcclm_[ln][cn] == (*fccl[l])[c]) {
                MIdx wn = bc.GetMIdx(cn);
                MIdx im = Unpack(fcimm_[ln][cn]);
                for (size_t d = 0; d < dim; ++d) { // periodic
                  (wn[d] < 0) && (im[d] += 1);
                  (wn[d] >= gs[d]) && (im[d] -= 1);
                }
                fcim_[l][c] = Pack(im);
                fndn = true;
              }
              if (fndn) break;
            }
            if (fndn) break;
          }
          if (!fndn) { // no same color, clear image
            fcim_[l][c] = Pack(MIdx(0));
          }
        }
      }
      for (auto lm : layers) { // check if old color disappeared
        if (fcclm_[lm][c] == kClNone) continue;
        bool fnd = false;
        for (auto l : layers) {
          if (fcclm_[lm][c] == (*fccl[l])[c]) {
            fnd = true;
            break;
          }
        }
        if (!fnd) { // color disappeared
          fcim_[lm][c] = Pack(MIdx(0)); // clear image
        }
      }
    }
    for (auto l : layers) {
      m.Comm(&fcim_[l]);
    }
  }
  if (sem("rotate")) {
    for (auto l : layers) {
      fcimm_[l] = fcim_[l];
      fcclm_[l] = *fccl[l];
    }
  }
}

} // namespace solver
