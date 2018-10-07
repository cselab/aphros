#pragma once

#include <limits>
#include <cmath>

#include "solver/solver.h"

namespace solver {

// Spherical averages.
template <class M_>
class Sphavg {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  static constexpr size_t dim = M::dim;

 public:
  struct Sph {
    Vect x;  // center
    Scal r;  // radius
    Scal h;  // kernel width
  };
  struct Avg {
    Vect x;  // average of position
    Scal b;  // total weight
    Avg() : x(0.), b(0.) {} 
  };
  // Constructor.
  // edim: effective dimension, 2 or 3
  Sphavg(M& m, size_t edim) : m(m), edim_(edim) {}
  // Computes averages over spheres.
  // fcu: volume fraction [a]
  // fcv: velocity [a]
  // fcvm: velocity from previous time step [a]
  // dt: time between time steps
  void Update(const FieldCell<Scal>& fcu,
              const FieldCell<Vect>& fcv, const FieldCell<Vect>& fcvm, Scal dt,
              const std::vector<Sph>& ss);
  // Returns spheres from last Update()
  const std::vector<Sph>& GetSph() const {
    return ss_;
  }
  // Returns averages from last Update()
  const std::vector<Avg>& GetAvg() const {
    return aa_;
  }

 private:
  M& m;
  size_t edim_; // effective dimension, 2 or 3
  std::vector<Sph> ss_; // spheres from last Update()
  std::vector<Avg> aa_; // averages from last Update()

  // Triangular kernel with support [-1,1]
  static Scal Kernel0(const Scal r) {
    return std::max(0., 1. - std::abs(r));
  }
  // Returns kernel value.
  // x: target 
  // s: sphere
  static Scal Kernel(const Vect& x, const Sph& s) {
    return Kernel0((x.dist(s.x) - s.r) / s.h);
  }
  // Returns projection to standard domain [0, m.GetGlobalLength()]
  Vect GetStd(const Vect& x) const {
    Vect r = x;
    auto l = m.GetGlobalLength();
    for (size_t d = 0; d < dim; ++d) {
      r[d] = (r[d] > 0. ? std::fmod(r[d], l[d]) : -std::fmod(-r[d], l[d]));
    }
    return r;
  }
  // Returns indices of bounding box (inclusive)
  // with center in standard domain [0, m.GetGlobalSize())
  Rect<MIdx> GetBox(const Sph& s) const {
    auto h = m.GetCellSize();
    MIdx wr((s.r + s.h) / h + 1);
    MIdx wx((GetStd(s.x) + Vect(0.5)) / h);
    return Rect<MIdx>(wx - wr, wx + wr);
  }
  // Returns indices of bounding box (inclusive) of mesh m
  // with center in standard domain [0, m.GetGlobalSize())
  Rect<MIdx> GetBox() const {
    auto& bc = m.GetInBlockCells();
    return Rect<MIdx>(bc.GetBegin(), bc.GetEnd() - MIdx(1));
  }
  // Returns 1 if segments [b,e] and [bo,eo] intersect
  static bool Inter(IntIdx b, IntIdx e, IntIdx bo, IntIdx eo) {
    return std::max(b, bo) <= std::min(e, eo);
  }
  // Returns 1 if rectangles r and ro intersect.
  // r: must lie in standard domain
  // ro: periodic conditions applied
  bool Inter(const Rect<MIdx>& r, const Rect<MIdx>& ro) const {
    MIdx b = r.lb;
    MIdx e = r.rt;
    MIdx bo = ro.lb;
    MIdx eo = ro.rt;
    auto gs = m.GetGlobalSize();
    // equivalent to intersection of projections in all directions
    bool p = true;
    for (size_t d = 0; d < dim; ++d) {
      bool q = false;
      q = q || Inter(b[d], bo[d], e[d], eo[d]);
      q = q || Inter(b[d], bo[d] + gs[d], e[d], eo[d] + gs[d]);
      q = q || Inter(b[d], bo[d] - gs[d], e[d], eo[d] - gs[d]);
      p = p && q;
    }
    return p;
  }
  void ClearAvg(size_t n) {
    aa_.clear();
    aa_.resize(n);
  }
};

template <class M_>
void Sphavg<M_>::Update(
    const FieldCell<Scal>& fcu,
    const FieldCell<Vect>& fcv, const FieldCell<Vect>& fcvm, Scal dt,
    const std::vector<Sph>& ss) {
  auto sem = m.GetSem("upd");

  auto& bc = m.GetIndexCells();

  if (sem("calc")) {
    ss_ = ss;
    ClearAvg(ss_.size());

    auto rm = GetBox();

    for (auto& s : ss) {
      auto rs = GetBox(s);
      if (Inter(rm, rs)) {
        std::cout << "inter: " 
            << "[" << rm.lb << "," << rm.rt << "]"
            << " "
            << "[" << rs.lb << "," << rs.rt << "]"
            ;
      }
    }
  }
  if (sem("reduce")) {
  }
}


} // namespace solver
