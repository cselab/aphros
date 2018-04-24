#pragma once

#include <exception>
#include <fstream>
#include <array>
#include <memory>

#include "mesh.hpp"
#include "linear.hpp"
#include "solver.hpp"
#include "advection.hpp"
#include "block.h"

namespace solver {

using namespace geom;

template <class Scal>
inline void Clip(Scal& a, Scal l, Scal u) {
  a = std::max(l, std::min(u, a));
}

template <class Scal>
inline void Clip(Scal& a) {
  Clip(a, 0., 1.);
};

// GetLineA() helper
// assuming 0 < u < 0.5, 0 < nx < ny
template <class Scal>
inline Scal GetLineA0(Scal nx, Scal ny, Scal u) {
  Scal u1 = 0.5 * nx / ny;
  if (u <= u1) {
    return -0.5 * (nx + ny) + std::sqrt(2. * nx * ny * u);
  } else {
    return ny * (u - 0.5);
  }
}

// Line constant by volume fraction in unit cell.
// n : normal
// u: volume fraction
// Returns:
// a: line constant
// Equation of reconstructed line 
// x.dot(n) = a
template <class Scal>
inline Scal GetLineA(const GVect<Scal, 3>& n, Scal u) {
  using Vect = GVect<Scal, 3>;

  Scal nx = std::abs(n[0]);
  Scal ny = std::abs(n[1]);
  if (ny < nx) {
    std::swap(nx, ny);
  }
  
  Clip(u);

  if (u < 0.5) {
    return GetLineA0(nx, ny, u);
  } else {
    return -GetLineA0(nx, ny, 1. - u);
  }
}

// Line constant by volume fraction in rectangular cell.
// n : normal
// u: volume fraction
// h: cell size
// Returns:
// a: line constant
// Equation of reconstructed line 
// x.dot(n) = a
template <class Scal>
inline Scal GetLineA(const GVect<Scal, 3>& n, Scal u, 
                     const GVect<Scal, 3>& h) {
  // TODO: check that GetLineA() is homogeneous wrt n
  return GetLineA(n * h, u);
}

// GetLineU() helper
// assuming a < 0, 0 < nx < ny
template <class Scal>
inline Scal GetLineU0(Scal nx, Scal ny, Scal a) {
  Scal a1 = 0.5 * (nx - ny);
  if (a <= a1) {
    if (nx == 0.) { // TODO revise
      return 0.;
    } else {
      return std::pow(a + 0.5 * (nx + ny), 2) / (2. * nx * ny);
    }
  } else {
    return 0.5 + a / ny;
  }
}

// Volume fraction from line constant in unit cell.
// n : normal
// a: line constant
// Returns:
// u: volume fraction
// Equation of reconstructed line 
// x.dot(n) = a
template <class Scal>
inline Scal GetLineU(const GVect<Scal, 3>& n, Scal a) {
  using Vect = GVect<Scal, 3>;

  Scal nx = std::abs(n[0]);
  Scal ny = std::abs(n[1]);
  if (ny < nx) {
    std::swap(nx, ny);
  }
  
  Clip(a, -0.5 * (nx + ny), 0.5 * (nx + ny));

  if (a < 0.) {
    return GetLineU0(nx, ny, a);
  } else {
    return 1. - GetLineU0(nx, ny, -a);
  }
}

// Volume fraction from line constant in rectangular cell.
// n : normal
// a: line constant
// h: cell size
// Returns:
// u: volume fraction
template <class Scal>
inline Scal GetLineU(const GVect<Scal, 3>& n, Scal a, 
                     const GVect<Scal, 3>& h) {
  // TODO: check that GetLineU() is homogeneous wrt n,a
  return GetLineU(n * h, a);
}

// GetLineVolX() helper
// assume dx > 0
template <class Scal>
inline Scal GetLineVolX0(const GVect<Scal, 3>& n, Scal a, 
                        const GVect<Scal, 3>& h, Scal dx) {
  using Vect = GVect<Scal, 3>;
  // Acceptor is a rectangular box adjacent to current cell.
  Vect hh(dx, h[1], h[2]); // acceptor size
  // Line constant for line advected by dx 
  // with origin at acceptor center
  // (e.g. shift 0 if dx=h[0], shift h[0]*0.5 if dx=0)
  Vect dc = Vect((h[0] - dx) *  0.5, 0., 0.); // shift of center
  Scal aa = a - n.dot(dc); // new line constant
  Scal uu = solver::GetLineU(n, aa, hh); // volume fraction
  Scal vv = hh.prod(); // acceptor volume
  Scal r = uu * vv;  // result
  r = std::min(r, GetLineU(n, a, h) * h.prod()); // limit by fluid in cell
  return r;
}

// Volume surplus in downwind adjacent cell after advection in x.
// (right if dx > 0, left if dx < 0)
// n : normal
// a: line constant
// h: cell size
// dx: advection displacement in x
// Returns:
// volume surplus in adjacent cell
template <class Scal>
inline Scal GetLineVolX(GVect<Scal, 3> n, Scal a, 
                        const GVect<Scal, 3>& h, Scal dx) {
  if (dx > 0.) {
    return GetLineVolX0(n, a, h, dx);
  } else {
    n[0] = -n[0];
    dx = -dx;
    return GetLineVolX0(n, a, h, dx);
  }
}

// Same as GetLineVolX in y.
template <class Scal>
inline Scal GetLineVolY(GVect<Scal, 3> n, Scal a, 
                        GVect<Scal, 3> h, Scal dy) {
  std::swap(n[0], n[1]);
  std::swap(h[0], h[1]);
  return GetLineVolX(n, a, h, dy);
}

// Fluid volume flux to downwind adjacent cell in x.
// n : normal
// h: cell size
// vx: mixture volume flux
// dt: time step
// Returns:
// fluid volume flux
template <class Scal>
inline Scal GetLineFluxX(const GVect<Scal, 3>& n, Scal a, 
                         const GVect<Scal, 3>& h, Scal vx, Scal dt) {
  Scal s = h[1] * h[2];  // face area
  Scal dx = vx / s * dt; // displacement
  Scal v = GetLineVolX(n, a, h, dx);
  if (vx < 0.) {
    v = -v;
  }
  return v / dt;
}

// Same as GetLineFluxX in y.
template <class Scal>
inline Scal GetLineFluxY(const GVect<Scal, 3>& n, Scal a, 
                         const GVect<Scal, 3>& h, Scal vy, Scal dt) {
  Scal s = h[0] * h[2];  // face area
  Scal dy = vy / s * dt; // displacement
  Scal v = GetLineVolY(n, a, h, dy);
  if (vy < 0.) {
    v = -v;
  }
  return v / dt;
}


template <class M>
class Vof : public AdvectionSolver<M> {
  using Mesh = M;
  using P = AdvectionSolver<Mesh>;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  static constexpr size_t dim = Mesh::dim;

  using P::m;
  using P::ffv_;
  using P::fcs_;
  LayersData<FieldCell<Scal>> fc_u_;
  MapFace<std::shared_ptr<ConditionFace>> mf_u_cond_;

  FieldCell<Scal> fc_a_; // alpha (plane constant)
  FieldCell<Vect> fc_n_; // n (normal to plane)
  FieldCell<Scal> fc_us_; // smooth field
  FieldFace<Scal> ff_fu_; // volume flux

 public:
  struct Par {
  };
  std::shared_ptr<Par> par;
  Par* GetPar() { return par.get(); }
  Vof(
      Mesh& m,
      const FieldCell<Scal>& fc_u_initial,
      const MapFace<std::shared_ptr<ConditionFace>>& mf_u_cond_,
      const FieldFace<Scal>* ff_volume_flux,
      const FieldCell<Scal>* fc_source,
      double t, double dt, std::shared_ptr<Par> par)
      : AdvectionSolver<Mesh>(t, dt, m, ff_volume_flux, fc_source)
      , mf_u_cond_(mf_u_cond_)
      , par(par)
      , fc_a_(m, 0)
      , fc_n_(m, Vect(0))
      , fc_us_(m, 0)
      , ff_fu_(m, 0) 
  {
    fc_u_.time_curr = fc_u_initial;
    Reconst(fc_u_.time_curr);
  }
  void StartStep() override {
    this->ClearIter();
    fc_u_.time_prev = fc_u_.time_curr;
    fc_u_.iter_curr = fc_u_.time_prev;
  }
  Scal Maxmod(Scal a, Scal b) {
    return std::abs(b) < std::abs(a) ? a : b;
  }
  void CalcNormal(const FieldCell<Scal>& uc) {
    auto uf = Interpolate(uc, mf_u_cond_, m);
    auto gc = Gradient(uf, m);
    for (auto c : m.AllCells()) {
      Vect g = gc[c];
      g[0] = Maxmod(g[0], 1e-10); // TODO: revise GetLine* to avoid this
      g[1] = Maxmod(g[1], 1e-10);
      fc_n_[c] = g / (-g.norm());
    }
  }
  // 2D
  void CalcNormalHeight(const FieldCell<Scal>& uc) {
    using MIdx = typename Mesh::MIdx;
    using Dir = typename Mesh::Dir;
    auto& bc = m.GetBlockCells();
    const int sw = 1; // stencil width
    const int sn = sw * 2 + 1; // stencil size
    GBlock<IdxCell, dim> bo(MIdx(-sw, -sw, 0), MIdx(sn, sn, 1)); // offset
    std::array<Scal, sn> e; // h[e]ight function
    auto h = GetCellSize();
    for (auto c : m.SuCells()) {
      Vect tn; // bes[t] normal
      Scal tk; // bes[t] slope with minimal abs
      for (Dir d : {Dir::i, Dir::j}) {
        // zero e
        for (size_t i = 0; i < sn; ++i) {
          e[i] = 0.;
        }
        auto w = bc.GetMIdx(c);
        // calc e
        for (auto o : bo) {
          e[o[size_t(d)] + sw] += uc[bc.GetIdx(w + o)];
        }
        // slope
        Scal km = (e[sw] - e[0]); // backward (minus)
        Scal kc = (e[sw + 1] - e[0]) / 2.; // centered
        Scal kp = (e[sw + 1] - e[sw]); // forward (plus)
        // best slope with maximum abs
        Scal k = Maxmod(km, Maxmod(kc, kp));
        // direction perpendicular 
        Dir dp(1 - size_t(d)); 
        // sign in dp
        Scal sg = uc[bc.GetIdx(w + MIdx(dp))] - uc[bc.GetIdx(w - MIdx(dp))];
        // normal
        Vect n;
        n[size_t(d)] = -k;
        n[size_t(dp)] = sg > 0. ? -1. : 1.;
        // check best
        if (d == Dir::i || (d == Dir::j && std::abs(k) < std::abs(tk))) {
          tn = n;
          tk = k;
        } 
      }
      tn[0] = Maxmod(tn[0], 1e-12); // TODO: revise GetLine* to avoid this
      tn[1] = Maxmod(tn[1], 1e-12);
      fc_n_[c] = tn / tn.norm();
    }
  }
  void Reconst(const FieldCell<Scal>& uc) {
    //CalcNormal(uc);
    CalcNormalHeight(uc);
    auto h = GetCellSize();
    for (auto c : m.AllCells()) {
      fc_a_[c] = GetLineA(fc_n_[c], uc[c], h);
    }
  }
  void Print(const FieldFace<Scal>& ff, std::string name) {
    using MIdx = typename Mesh::MIdx;
    auto ibc = m.GetInBlockCells();
    auto bc = m.GetBlockCells();
    auto bf = m.GetBlockFaces();
    MIdx wb = ibc.GetBegin();
    MIdx we = ibc.GetEnd();
    MIdx wc = (wb + we - MIdx(1)) / 2;
    std::cerr << std::setw(5) << name << " = ";
    for (int i = wb[1]; i <= we[1]; ++i) {
      MIdx w = wc;
      w[1] = i;
      using Dir = typename Mesh::Dir;
      IdxFace f = bf.GetIdx(w, Dir::j);
      std::cerr << std::setw(10) << ff[f] << " ";
    }
    std::cerr << std::endl;
  }

  void Print(const FieldCell<Scal>& fc, std::string name) {
    using MIdx = typename Mesh::MIdx;
    auto ibc = m.GetInBlockCells();
    auto bc = m.GetBlockCells();
    MIdx wb = ibc.GetBegin();
    MIdx we = ibc.GetEnd();
    MIdx wc = (wb + we - MIdx(1)) / 2;
    std::cerr << std::setw(5) << name << " = ";
    for (int i = wb[1]; i < we[1]; ++i) {
      MIdx w = wc;
      w[1] = i;
      IdxCell c = bc.GetIdx(w);
      std::cerr << std::setw(10) << fc[c] << " ";
    }
    std::cerr << std::endl;
  }

  Vect GetCellSize() const {
    Vect h; // cell size
    // XXX: Adhoc for structured 3D mesh
    IdxCell c0(0);
    h = m.GetNode(m.GetNeighbourNode(c0, 7)) - 
        m.GetNode(m.GetNeighbourNode(c0, 0));
    assert(std::abs(h.prod() - m.GetVolume(c0)) < 1e-10);
    return h;
  }

  void MakeIteration() override {
    auto sem = m.GetSem("iter");
    if (sem("init")) {
      auto& uc = fc_u_.iter_curr;
      const Scal dt = this->GetTimeStep();
      for (auto c : m.Cells()) {
        uc[c] = fc_u_.time_prev[c] +  // previous time step
            dt * (*fcs_)[c]; // source
      }
    }
    const bool split = true;
    for (size_t d = 0; d < (split ? dim : 1); ++d) {
      auto& uc = fc_u_.iter_curr;
      if (sem("adv")) {
        auto h = GetCellSize();

        auto& ffv = *ffv_; // [f]ield [f]ace [v]olume flux
        const Scal dt = this->GetTimeStep();
        for (auto f : m.Faces()) {
          Vect vd(0); 
          vd[d] = 1.;

          Scal v = ffv[f] * vd.dot(m.GetNormal(f)); // mixture flux
          bool p = (v > 0.);
          IdxCell c = m.GetNeighbourCell(f, p ? 0 : 1); // upwind cell
          if (d == 0) {
            ff_fu_[f] = GetLineFluxX(fc_n_[c], fc_a_[c], h, v, dt);
          } else if (d == 1) {
            ff_fu_[f] = GetLineFluxY(fc_n_[c], fc_a_[c], h, v, dt);
          } else { // d == 3
            ff_fu_[f] = 0.;
          }
        }

        for (auto c : m.Cells()) {
          Scal s = 0.; // sum of fluxes
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetNeighbourFace(c, q);
            s += ff_fu_[f] * m.GetOutwardFactor(c, q);
          }

          uc[c] += dt / m.GetVolume(c) * (-s);
        }

        m.Comm(&uc);
      }
      if (sem("reconst")) {
        Reconst(uc);
      }
    }

    if (sem("stat")) {
      this->IncIter();
    }
  }
  void FinishStep() override {
    fc_u_.time_curr = fc_u_.iter_curr;
    this->IncTime();
  }
  const FieldCell<Scal>& GetField(Layers l) override {
    return fc_u_.Get(l);
  }
  const FieldCell<Scal>& GetAlpha() {
    return fc_a_;
  }
  const FieldCell<Vect>& GetNormal() {
    return fc_n_;
  }
  using P::GetField;
};

} // namespace solver
