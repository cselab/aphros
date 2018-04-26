#pragma once

#include <exception>
#include <fstream>
#include <array>
#include <memory>

#include "advection.h"
#include "geom/block.h"

namespace solver {

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

// GetLineA() helper for unit cell.
// n : normal
// u: volume fraction
// Returns:
// a: line constant
// Equation of reconstructed line 
// x.dot(n) = a
template <class Scal>
inline Scal GetLineA1(const GVect<Scal, 3>& n, Scal u) {
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
  return GetLineA1(n * h, u);
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

// GetLineU() helper for unit cell
// n : normal
// a: line constant
// Returns:
// u: volume fraction
// Equation of reconstructed line 
// x.dot(n) = a
template <class Scal>
inline Scal GetLineU1(const GVect<Scal, 3>& n, Scal a) {
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
  return GetLineU1(n * h, a);
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
  Scal uu = GetLineU(n, aa, hh); // volume fraction
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

// GetLineVolStrX() helper
// assume dx > 0
template <class Scal>
inline Scal GetLineVolStrX0(const GVect<Scal, 3>& n, Scal a, 
                            const GVect<Scal, 3>& h, Scal dx, Scal dxu) {
  using Vect = GVect<Scal, 3>;
  Scal u = GetLineU(n, a, h); // volume fraction
  Vect sh(h[0] + dx - dxu, h[1], h[2]); // stretched size
  Vect sn = n / sh; // stretched normal
  Scal sa = GetLineA(sn, u, sh); // stretched line constant
  Vect dc = Vect(dxu, 0., 0.); // shift of center
  return GetLineVolX0(sn, sa, sh, dx);
}

// Volume surplus in downwind adjacent cell after stretching in x
// (right if dx > 0, left if dx < 0)
// n : normal
// a: line constant
// h: cell size
// dx: displacement in x
// dxu: displacement in x of the other face upwind
// Returns:
// volume surplus in adjacent cell
template <class Scal>
inline Scal GetLineVolStrX(GVect<Scal, 3> n, Scal a, 
                        const GVect<Scal, 3>& h, Scal dx, Scal dxu) {
  if (dx > 0.) {
    return GetLineVolStrX0(n, a, h, dx, dxu);
  } else {
    n[0] = -n[0];
    dx = -dx;
    dxu = -dxu;
    return GetLineVolStrX0(n, a, h, dx, dxu);
  }
}

// Same as GetLineVolStrX in y.
template <class Scal>
inline Scal GetLineVolStrY(GVect<Scal, 3> n, Scal a, 
                           GVect<Scal, 3> h, Scal dy, Scal dyu) {
  std::swap(n[0], n[1]);
  std::swap(h[0], h[1]);
  return GetLineVolStrX(n, a, h, dy, dyu);
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

// Fluid volume flux to downwind adjacent cell in x after stretching.
// n : normal
// h: cell size
// vx: mixture volume flux
// vxu: mixture volume flux upwind face
// dt: time step
// Returns:
// fluid volume flux
template <class Scal>
inline Scal GetLineFluxStrX(const GVect<Scal, 3>& n, Scal a, 
                            const GVect<Scal, 3>& h, 
                            Scal vx, Scal vxu, Scal dt) {
  Scal s = h[1] * h[2];  // face area
  Scal dx = vx / s * dt; // displacement
  Scal dxu = vxu / s * dt; // displacement on upwind face
  Scal v = GetLineVolStrX(n, a, h, dx, dxu);
  if (vx < 0.) {
    v = -v;
  }
  return v / dt;
}

// Same as GetLineFluxStrX in y.
template <class Scal>
inline Scal GetLineFluxStrY(const GVect<Scal, 3>& n, Scal a, 
                            const GVect<Scal, 3>& h, 
                            Scal vy, Scal vyu, Scal dt) {
  Scal s = h[0] * h[2];  // face area
  Scal dy = vy / s * dt; // displacement
  Scal dyu = vyu / s * dt; // displacement
  Scal v = GetLineVolStrY(n, a, h, dy, dyu);
  if (vy < 0.) {
    v = -v;
  }
  return v / dt;
}

template <class M_>
class Vof : public AdvectionSolver<M_> {
  using M = M_;
  using P = AdvectionSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;

  using P::m;
  using P::ffv_;
  using P::fcs_;
  LayersData<FieldCell<Scal>> fc_u_;
  MapFace<std::shared_ptr<ConditionFace>> mf_u_cond_;

  FieldCell<Scal> fc_a_; // alpha (plane constant)
  FieldCell<Vect> fc_n_; // n (normal to plane)
  FieldCell<Scal> fc_us_; // smooth field
  FieldFace<Scal> ff_fu_; // volume flux
  size_t count_ = 0; // number of MakeIter() calls, used for splitting

 public:
  struct Par {
  };
  std::shared_ptr<Par> par;
  Par* GetPar() { return par.get(); }
  Vof(
      M& m,
      const FieldCell<Scal>& fc_u_initial,
      const MapFace<std::shared_ptr<ConditionFace>>& mf_u_cond_,
      const FieldFace<Scal>* ff_volume_flux,
      const FieldCell<Scal>* fc_source,
      double t, double dt, std::shared_ptr<Par> par)
      : AdvectionSolver<M>(t, dt, m, ff_volume_flux, fc_source)
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
      fc_n_[c] = g;
    }
  }
  // 2D
  void CalcNormalHeight(const FieldCell<Scal>& uc) {
    using MIdx = typename M::MIdx;
    using Dir = typename M::Dir;
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
        Scal kc = (e[sw + 1] - e[0]) * 0.5; // centered
        Scal kp = (e[sw + 1] - e[sw]); // forward (plus)
        // best slope with maximum abs
        //Scal k = Maxmod(km, Maxmod(kc, kp));
        Scal k = kc; // XXX: force centered approx
        // direction perpendicular 
        Dir dp(1 - size_t(d)); 
        // sign in dp
        Scal sg = uc[bc.GetIdx(w + MIdx(dp))] - uc[bc.GetIdx(w - MIdx(dp))];
        // normal
        Vect n;
        n[size_t(d)] = -k;
        n[size_t(dp)] = sg > 0. ? -1. : 1.;
        // check best with minimal abs
        if (d == Dir::i || (d == Dir::j && std::abs(k) < std::abs(tk))) {
          tn = n;
          tk = k;
        } 
      }
      fc_n_[c] = tn;
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
    using MIdx = typename M::MIdx;
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
      using Dir = typename M::Dir;
      IdxFace f = bf.GetIdx(w, Dir::j);
      std::cerr << std::setw(10) << ff[f] << " ";
    }
    std::cerr << std::endl;
  }

  void Print(const FieldCell<Scal>& fc, std::string name) {
    using MIdx = typename M::MIdx;
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
    using Dir = typename M::Dir;
    using MIdx = typename M::MIdx;
    auto& bf = m.GetBlockFaces();
    std::vector<Dir> dd;
    if (count_ % 2 == 0) {
      dd = {Dir::i, Dir::j};
    } else {
      dd = {Dir::j, Dir::i};
    }
    for (Dir d : dd) {
      auto& uc = fc_u_.iter_curr;
      auto& bc = m.GetBlockCells();
      auto& bf = m.GetBlockFaces();
      if (sem("adv")) {
        auto h = GetCellSize();
        auto& ffv = *ffv_; // [f]ield [f]ace [v]olume flux
        const Scal dt = this->GetTimeStep();
        for (auto c : m.Cells()) {
          auto w = bc.GetMIdx(c);
          const Scal lc = m.GetVolume(c);
          // faces
          IdxFace fm = bf.GetIdx(w, d);
          IdxFace fp = bf.GetIdx(w + MIdx(d), d);
          // mixture volume fluxes
          const Scal vm = ffv[fm];
          const Scal vp = ffv[fp];
          // mixture volume cfl
          const Scal sm = vm * dt / lc;
          const Scal sp = vp * dt / lc;
          const Scal ds = sp - sm;
          // upwind cells
          IdxCell cum = m.GetNeighbourCell(fm, vm > 0. ? 0 : 1);
          IdxCell cup = m.GetNeighbourCell(fp, vp > 0. ? 0 : 1);
          // phase 1 volume fluxes
          Scal qm, qp; 
          if (d == dd[0]) { // Euler Implicit
            if (d == Dir::i) {
              qm = GetLineFluxX(fc_n_[cum], fc_a_[cum], h, vm, dt);
              qp = GetLineFluxX(fc_n_[cup], fc_a_[cup], h, vp, dt);
            } else if (d == Dir::j) {
              qm = GetLineFluxY(fc_n_[cum], fc_a_[cum], h, vm, dt);
              qp = GetLineFluxY(fc_n_[cup], fc_a_[cup], h, vp, dt);
            } else if (d == Dir::k) {
              // nop
            }
            // phase 1 cfl
            const Scal lm = qm * dt / lc;
            const Scal lp = qp * dt / lc;
            const Scal dl = lp - lm;
            uc[c] = (uc[c] - dl) / (1. - ds);
          } else { // Lagrange Explicit
            // upwind faces
            IdxFace fum = bf.GetIdx(vm > 0. ? w - MIdx(d) : w, d);
            IdxFace fup = bf.GetIdx(vp > 0. ? w : w + MIdx(d), d);
            // upwind fluxes
            Scal vum = ffv[fum];
            Scal vup = ffv[fup];
            if (d == Dir::i) {
              qm = GetLineFluxStrX(fc_n_[cum], fc_a_[cum], h, vm, vum, dt);
              qp = GetLineFluxStrX(fc_n_[cup], fc_a_[cup], h, vp, vup, dt);
            } else if (d == Dir::j) {
              qm = GetLineFluxStrY(fc_n_[cum], fc_a_[cum], h, vm, vum, dt);
              qp = GetLineFluxStrY(fc_n_[cup], fc_a_[cup], h, vp, vup, dt);
            } else if (d == Dir::k) {
              // nop
            }
            // phase 1 cfl
            const Scal lm = qm * dt / lc;
            const Scal lp = qp * dt / lc;
            const Scal dl = lp - lm;
            uc[c] = uc[c] * (1. + ds) - dl;
          }
          Clip(uc[c]);
        }
        m.Comm(&uc);
      }
      if (sem("reconst")) {
        Reconst(uc);
      }
    }

    if (sem("stat")) {
      this->IncIter();
      ++count_;
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
