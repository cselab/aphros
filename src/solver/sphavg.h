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
    Vect x = Vect(0);  // center
    Scal r = 1.;  // radius
    Scal h = 1.;  // kernel width
    Sph() = default;
    Sph(const Vect& x, Scal r, Scal h) : x(x), r(r), h(h) {} 
  };
  struct Avg {
    Scal b = 0.;  // total weight
    Vect x = Vect(0);  // average of position
    Vect v = Vect(0);  // velocity
    Vect dvt = Vect(0);  // dv/dt
    Vect dvx = Vect(0);  // dv/dx * v
    Scal p = 0.;
    Avg() = default;
    static std::vector<std::string> GetNames() {
      std::vector<std::string> r;
      auto as = [&r](std::string n) {
        r.push_back(n);
      };
      auto av = [&as](std::string n) {
        as(n + "x");
        as(n + "y");
        as(n + "z");
      };
      as("b");
      av("");
      av("v");
      av("dvt");
      av("dvx");
      as("p");
      return r;
    }
    void App(Scal a, std::vector<Scal>& g) const {
      g.push_back(a);
    }
    void App(const Vect& a, std::vector<Scal>& g) const {
      App(a[0], g);
      App(a[1], g);
      App(a[2], g);
    }
    std::vector<Scal> Ser() const {
      std::vector<Scal> g;
      App(b, g);
      App(x, g);
      App(v, g);
      App(dvt, g);
      App(dvx, g);
      App(p, g);
      return g;
    }
    void Ext(Scal& a, const std::vector<Scal>& g, size_t& i) {
      a = g[i++];
    }
    void Ext(Vect& a, const std::vector<Scal>& g, size_t& i) {
      Ext(a[0], g, i);
      Ext(a[1], g, i);
      Ext(a[2], g, i);
    }
    void Des(const std::vector<Scal>& g) {
      size_t i = 0;
      Ext(b, g, i);
      Ext(x, g, i);
      Ext(v, g, i);
      Ext(dvt, g, i);
      Ext(dvx, g, i);
      Ext(p, g, i);
    }
  };
  // Constructor.
  // edim: effective dimension, 2 or 3
  Sphavg(M& m, size_t edim) : m(m), edim_(edim) {}
  // Computes averages over spheres.
  // fcu: volume fraction [a]
  // fcv: velocity [a]
  // fcvm: velocity from previous time step [a]
  // fcp: pressure [a]
  // dt: time between time steps
  void Update(const FieldCell<Scal>& fcu,
              const FieldCell<Vect>& fcv, const FieldCell<Vect>& fcvm, Scal dt,
              const FieldCell<Scal>& fcp,
              const std::vector<Sph>& ss);
  // Returns spheres from last Update()
  const std::vector<Sph>& GetSph() const {
    return ss_;
  }
  // Returns averages from last Update()
  const std::vector<Avg>& GetAvg() const {
    return aa_;
  }
  const std::vector<std::string> GetNames() const {
    return Avg::GetNames();
  }

 public:
  M& m;
  size_t edim_; // effective dimension, 2 or 3
  std::vector<Sph> ss_; // spheres from last Update()
  std::vector<Avg> aa_; // averages from last Update()
  std::vector<std::vector<Scal>> vv_; // tmp for reduction

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
      r[d] = (r[d] >= 0. ? std::fmod(r[d], l[d]) : -std::fmod(-r[d], l[d]));
    }
    return r;
  }
  // Returns projection to standard domain [0, m.GetGlobalSize())
  MIdx GetStd(const MIdx& w) const {
    MIdx r = w;
    auto gs = m.GetGlobalSize();
    for (size_t d = 0; d < dim; ++d) {
      r[d] = (std::abs(r[d]) * gs[d] + r[d]) % gs[d];
    }
    return r;
  }
  // Returns indices of bounding box (inclusive)
  // with center in standard domain [0, m.GetGlobalSize())
  Rect<MIdx> GetBox(const Sph& s) const {
    auto h = m.GetCellSize();
    MIdx wr(Vect(s.r + s.h) / h + Vect(1.));
    for (size_t d = 0; d < dim; ++d) {
      wr[d] = std::min(wr[d], m.GetGlobalSize()[d]);
    }
    MIdx wx(GetStd(s.x) / h + Vect(0.5));
    if (edim_ < 3) {
      wr[2] = 0;
      wx[2] = 0;
    }
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
    for (size_t d = 0; d < dim; ++d) {
      bool q = false;
      q = q || Inter(b[d], e[d], bo[d], eo[d]);
      q = q || Inter(b[d], e[d], bo[d] + gs[d], eo[d] + gs[d]);
      q = q || Inter(b[d], e[d], bo[d] - gs[d], eo[d] - gs[d]);
      if (!q) {
        return false;
      }
    }
    return true;
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
    const FieldCell<Scal>& fcp,
    const std::vector<Sph>& ss) {
  auto sem = m.GetSem("upd");

  auto& bd = m.GetIndexCells();
  auto& bc = m.GetInBlockCells();
  (void) fcu;

  if (sem("calc")) {
    ss_ = ss;
    ClearAvg(ss_.size());
    auto h = m.GetCellSize();

    auto rm = GetBox();

    for (size_t i = 0; i < ss_.size(); ++i) {
      auto& s = ss_[i];
      auto& a = aa_[i];
      auto rs = GetBox(s);
      if (Inter(rm, rs)) {
        // traverse bounding box
        typename M::BlockCells bs(rs.lb, rs.GetDimensions() + MIdx(1));
        for (auto w : bs) {
          MIdx wt = GetStd(w);
          Sph st = s;
          st.x = GetStd(st.x);
          if (bc.IsInside(wt)) {
            Vect x(Vect(w) * h);
            Scal b = Kernel(x, st);
            IdxCell c = bd.GetIdx(wt);
            a.b += b;
            a.x += x * b;
            a.v += fcv[c] * b;
            a.dvt += (fcv[c] - fcvm[c]) * (b / dt);
            a.dvx += Vect(0) * b; // TODO: fill dvx
            a.p += fcp[c] * b;
          }
        }
      }
    }

    vv_.clear();
    // serialize
    for (size_t i = 0; i < ss_.size(); ++i) {
      vv_.push_back(aa_[i].Ser());
    }

    using TVS = typename M::template OpCatVT<Scal>; 
    m.Reduce(std::make_shared<TVS>(&vv_));
  }
  if (sem("reduce")) {
    // root has concatenation of all vv_
    if (m.IsRoot()) {
      if (vv_.size() % ss_.size()) {
        throw std::runtime_error(
            "reduce: vv_.size()=" + std::to_string(vv_.size()) +
            " notdiv ss_.size()=" + std::to_string(ss_.size()));
      }
      auto e = ss_.size();
      for (size_t i = 0; i < e; ++i) {
        for (size_t ii = i + e; ii < vv_.size(); ii += e) {
          for (size_t j = 0; j < vv_[i].size(); ++j) {
            vv_[i][j] += vv_[ii][j];
          }
        }
      }
      vv_.resize(e);
      for (size_t i = 0; i < e; ++i) {
        for (size_t j = 1; j < vv_[i].size(); ++j) {
          vv_[i][j] /= vv_[i][0];
        }
      }
    }
    using TVS = typename M::template OpCatVT<Scal>; 
    m.Bcast(std::make_shared<TVS>(&vv_));
  }
  if (sem("bcast")) {
    if (vv_.size() != ss_.size()) {
      throw std::runtime_error(
          "bcast: vv_.size()=" + std::to_string(vv_.size()) +
          " != ss_.size()=" + std::to_string(ss_.size()));
    }
    // deserialize
    for (size_t i = 0; i < ss_.size(); ++i) {
      aa_[i].Des(vv_[i]);
    }
  }
}

} // namespace solver
