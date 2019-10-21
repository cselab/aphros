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
      : m(m), layers(layers), fcim_(layers) {
    fcim_.InitAll(FieldCell<Scal>(m, Pack(MIdx(0))));
  }
  // fccl: current color
  // fcclm: previous color
  void Update(const Multi<const FieldCell<Scal>*>& fccl,
              const Multi<const FieldCell<Scal>*>& fcclm);
  // Returns image vector, number of passes through periodic boundaries
  MIdx GetImage(size_t l, IdxCell c) const {
    return Unpack(fcim_[l][c]);
  }

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
  Multi<FieldCell<Scal>> fcim_;  // image current [a]
};

template <class M_>
void Trackerm<M_>::Update(const Multi<const FieldCell<Scal>*>& fccl,
                          const Multi<const FieldCell<Scal>*>& fcclm) {
  auto sem = m.GetSem("trackerm");
  struct {
    Multi<FieldCell<Scal>> fcimm; // image previous [a]
  }* ctx(sem);
  auto& fcimm = ctx->fcimm;
  if (sem("update")) {
    fcimm = fcim_; // save previous
    MIdx gs = m.GetGlobalSize();
    auto& bc = m.GetIndexCells();

    for (auto c : m.Cells()) {
      for (auto l : layers) { // check if new color appeared
        if ((*fccl[l])[c] == kClNone) { // no color, clear image
          fcim_[l][c] = Pack(MIdx(0));
          continue;
        }
        bool fndm = false;
        for (auto lm : layers) {
          if ((*fcclm[lm])[c] == (*fccl[l])[c]) {
            fndm = true;
            break;
          }
        }
        if (!fndm) { // new color, find same color in neighbors
          bool fndn = false;
          for (auto q : m.Nci(c)) {
            auto cn = m.GetCell(c, q);
            for (auto ln : layers) {
              if ((*fcclm[ln])[cn] == (*fccl[l])[c]) {
                MIdx wn = bc.GetMIdx(cn);
                MIdx im = Unpack(fcimm[ln][cn]);
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
          if (!fndn) { // no same color in neighbors, clear image
            fcim_[l][c] = Pack(MIdx(0));
          }
        }
      }
      for (auto lm : layers) { // check if old color disappeared
        if ((*fcclm[lm])[c] == kClNone) continue;
        bool fnd = false;
        for (auto l : layers) {
          if ((*fcclm[lm])[c] == (*fccl[l])[c]) {
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
}

} // namespace solver
