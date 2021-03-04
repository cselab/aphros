// Created by Petr Karnakov on 07.10.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cmath>
#include <limits>

#include "solver/solver.h"

// Spherical averages.
template <class M_>
class Sphavg {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Vect3 = generic::Vect<Scal, 3>;
  using MIdx = typename M::MIdx;
  static constexpr size_t dim = M::dim;

 public:
  struct Sph {
    Vect x = Vect(0); // center
    Scal r = 1; // radius
    Scal h = 1; // kernel width
    Sph() = default;
    Sph(const Vect& x_, Scal r_, Scal h_) : x(x_), r(r_), h(h_) {}
  };
  struct Avg {
    Scal b = 0; // total weight
    Vect x = Vect(0); // average of position
    Vect v = Vect(0); // velocity
    Vect dvt = Vect(0); // dv/dt
    Vect dvx = Vect(0); // dv/dx * v
    Vect dvmat = Vect(0); // dv/dt + dv/dx * v
    Vect gp = Vect(0); // gp/dx
    Scal vm = 0; // velocity magnitude
    Scal p = 0; // pressure
    Scal r = 1; // equivalent radius
    Scal rhm = 1; // shell inner
    Scal rhp = 1; // shell outer
    Vect gvx = Vect(0); // dv/dx
    Vect gvy = Vect(0); // dv/dy
    Vect gvz = Vect(0); // dv/dz
    Avg() = default;
    static std::vector<std::string> GetNames() {
      std::vector<std::string> r;
      auto as = [&r](std::string n) { r.push_back(n); };
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
      av("dvmat");
      av("gp");
      as("vm");
      as("p");
      av("gvx");
      av("gvy");
      av("gvz");
      as("r");
      as("rhm");
      as("rhp");
      return r;
    }
    void Serialize(Scal a, std::vector<Scal>& g) const {
      g.push_back(a);
    }
    void Serialize(const Vect& a, std::vector<Scal>& g) const {
      for (auto d : M::dirs) {
        Serialize(a[d], g);
      }
    }
    std::vector<Scal> Serialize() const {
      std::vector<Scal> g;
      Serialize(b, g);
      Serialize(x, g);
      Serialize(v, g);
      Serialize(dvt, g);
      Serialize(dvx, g);
      Serialize(dvmat, g);
      Serialize(gp, g);
      Serialize(vm, g);
      Serialize(p, g);
      Serialize(gvx, g);
      Serialize(gvy, g);
      Serialize(gvz, g);
      return g;
    }
    std::vector<Scal> SerOut() const {
      std::vector<Scal> g;
      Serialize(b, g);
      Serialize(x, g);
      Serialize(v, g);
      Serialize(dvt, g);
      Serialize(dvx, g);
      Serialize(dvmat, g);
      Serialize(gp, g);
      Serialize(vm, g);
      Serialize(p, g);
      Serialize(gvx, g);
      Serialize(gvy, g);
      Serialize(gvz, g);
      Serialize(r, g);
      Serialize(rhm, g);
      Serialize(rhp, g);
      return g;
    }
    void Ext(Scal& a, const std::vector<Scal>& g, size_t& i) {
      a = g[i++];
    }
    void Ext(Vect& a, const std::vector<Scal>& g, size_t& i) {
      for (auto d : M::dirs) {
        Ext(a[d], g, i);
      }
    }
    void Deserialize(const std::vector<Scal>& g) {
      size_t i = 0;
      Ext(b, g, i);
      Ext(x, g, i);
      Ext(v, g, i);
      Ext(dvt, g, i);
      Ext(dvx, g, i);
      Ext(dvmat, g, i);
      Ext(gp, g, i);
      Ext(vm, g, i);
      Ext(p, g, i);
      Ext(gvx, g, i);
      Ext(gvy, g, i);
      Ext(gvz, g, i);
    }
  };
  // Constructor.
  // edim: effective dimension, 2 or 3
  Sphavg(M& m_, size_t edim) : m(m_), edim_(edim) {}
  // Computes averages over spheres.
  // fcu: volume fraction [a]
  // fcv: velocity [a]
  // fcvm: velocity from previous time step [a]
  // fcp: pressure [a]
  // dt: time between time steps
  void Update(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcv,
      const FieldCell<Vect>& fcvm, Scal dt, const FieldCell<Scal>& fcp,
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
    if (edim_ < 3 && M::dim > 2) {
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
    MIdx b = r.low;
    MIdx e = r.high;
    MIdx bo = ro.low;
    MIdx eo = ro.high;
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
    const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcv,
    const FieldCell<Vect>& fcvm, Scal dt, const FieldCell<Scal>& fcp,
    const std::vector<Sph>& ss) {
  auto sem = m.GetSem("upd");

  auto& bd = m.GetIndexCells();
  auto& bc = m.GetInBlockCells();
  (void)fcu;

  if (sem("calc")) {
    ss_ = ss;
    ClearAvg(ss_.size());
    auto h = m.GetCellSize();

    auto rm = GetBox();

    // derivative in direction d at point w
    auto sd = [&h, &bd](const FieldCell<Scal>& f, MIdx w, size_t d) {
      MIdx wp = w;
      ++wp[d];
      MIdx wm = w;
      --wm[d];
      IdxCell cp = bd.GetIdx(wp);
      IdxCell cm = bd.GetIdx(wm);
      return (f[cp] - f[cm]) / (2. * h[d]);
    };
    // derivative in direction d at point w
    auto sdv = [&h, &bd](const FieldCell<Vect>& f, MIdx w, size_t d) {
      MIdx wp = w;
      ++wp[d];
      MIdx wm = w;
      --wm[d];
      IdxCell cp = bd.GetIdx(wp);
      IdxCell cm = bd.GetIdx(wm);
      return (f[cp] - f[cm]) / (2. * h[d]);
    };
    /*
    // derivative of component p in direction d at point w
    auto vd = [&h,&bd](const FieldCell<Vect>& f, size_t p, MIdx w, size_t d) {
      MIdx wp = w;
      ++wp[d];
      MIdx wm = w;
      --wm[d];
      IdxCell cp = bd.GetIdx(wp);
      IdxCell cm = bd.GetIdx(wm);
      return (f[cp][p] - f[cm][p]) / (2. * h[d]);
    };
    */

    for (size_t i = 0; i < ss_.size(); ++i) {
      auto& s = ss_[i];
      auto& a = aa_[i];
      auto rs = GetBox(s);
      if (Inter(rm, rs)) {
        // traverse bounding box
        typename M::BlockCells bs(rs.low, rs.GetDimensions() + MIdx(1));
        for (auto w : bs) {
          MIdx wt = GetStd(w);
          Sph st = s;
          Vect xd = st.x - GetStd(st.x); // shift from standard
          st.x = GetStd(st.x);
          if (bc.IsInside(wt)) {
            Vect x(Vect(w) * h);
            Scal b = Kernel(x, st); // weight
            IdxCell c = bd.GetIdx(wt);
            auto v = fcv[c]; // velocity
            a.b += b;
            a.x += (x + xd) * b;
            a.v += v * b;
            a.dvt += (v - fcvm[c]) * (b / dt);

            for (size_t d = 0; d < dim; ++d) {
              a.dvx += sdv(fcv, wt, d) * (v[d] * b);
            }
            for (size_t p = 0; p < dim; ++p) {
              a.gp[p] += sd(fcp, wt, p) * b;
            }
            a.dvmat = a.dvt + a.dvx;
            a.p += fcp[c] * b;
            a.vm += fcv[c].norm() * b;
            // velocity gradient
            {
              auto tx = sdv(fcv, wt, 0) * b;
              auto ty = sdv(fcv, wt, 1) * b;
              auto tz = sdv(fcv, wt, 2) * b;
              a.gvx += Vect(Vect3(tx[0], ty[0], tz[0]));
              a.gvy += Vect(Vect3(tx[1], ty[1], tz[1]));
              if (M::dim > 2) {
                a.gvz += Vect(Vect3(tx[2], ty[2], tz[2]));
              }
            }
          }
        }
      }
    }

    vv_.clear();
    // serialize
    for (size_t i = 0; i < ss_.size(); ++i) {
      vv_.push_back(aa_[i].Serialize());
    }

    m.Reduce(&vv_, Reduction::concat);
  }
  if (sem("reduce")) {
    // root has concatenation of all vv_
    if (m.IsRoot()) {
      fassert(
          vv_.size() % ss_.size() == 0,
          "reduce: vv_.size()=" + std::to_string(vv_.size()) +
              " not divisible by ss_.size()=" + std::to_string(ss_.size()));
      auto e = ss_.size();
      for (size_t i = 0; i < e; ++i) {
        // append all to [i]
        for (size_t ii = i + e; ii < vv_.size(); ii += e) {
          for (size_t j = 0; j < vv_[i].size(); ++j) {
            vv_[i][j] += vv_[ii][j];
          }
        }
      }
      vv_.resize(e);
      for (size_t i = 0; i < e; ++i) {
        for (size_t j = 1; j < vv_[i].size(); ++j) {
          vv_[i][j] /= vv_[i][0]; // assume weight is [0]
        }
      }
    }
    m.Bcast(&vv_);
  }
  if (sem("bcast")) {
    fassert_equal(vv_.size(), ss_.size());
    // deserialize
    for (size_t i = 0; i < ss_.size(); ++i) {
      aa_[i].Deserialize(vv_[i]);
    }
    // append r
    for (size_t i = 0; i < ss_.size(); ++i) {
      auto& a = aa_[i];
      auto& s = ss_[i];
      a.r = s.r;
      a.rhm = s.r - s.h;
      a.rhp = s.r + s.h;
    }
  }
}
