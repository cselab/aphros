#pragma once

#include <exception>
#include <fstream>
#include <array>
#include <memory>

#include "advection.h"
#include "geom/block.h"
#include "dump/dumper.h"

namespace solver {

template <class Scal>
inline void Clip(Scal& a, Scal l, Scal u) {
  a = std::max(l, std::min(u, a));
}

template <class Scal>
inline void Clip(Scal& a) {
  Clip(a, 0., 1.);
}

template <class T>
inline T cube(T a) {
  return a * a * a;
}

template <class T>
inline T sqr(T a) {
  return a * a;
}

template <class Scal>
Scal SolveCubic(Scal a, Scal b, Scal c, Scal d) {
  Scal p = (3. * a * c - b * b) / (3. * a * a);
  Scal q = (2. * cube(b) - 9. * a * b * c + 27. * a * a * d) / (27. * cube(a));
  int k = 1; // 0,1,2
  Scal t1 = 2. * std::sqrt(-p / 3.) *
      std::cos(
          1. / 3. * 
          std::acos(3. * q * std::sqrt(-3. / p) / (2. * p)) - 
          2. * M_PI * k / 3.);
}

// GetLineA() helper
// assuming 0 < u < 0.5, 0 < nx < ny < nz
template <class Scal>
inline Scal GetLineA0(Scal nx, Scal ny, Scal nz, Scal u) {
  Scal f;

  if (6. * ny * nz * u <= sqr(nx)) {
    f = std::pow(6. * nx * ny * nz * u, 1. / 3.);
  } else if (6. * ny * nz * u <= 3. * sqr(ny) - 3 * nx * ny + sqr(nx)) {
    f = 0.5 * nx + std::sqrt(2. * ny * nz * u - sqr(nx) / 12.);
  } else if (2. * nz * u < nx + ny && 
      6. * nx * ny * nz * u < -cube(nz) + 3. * sqr(nz) * (nx + ny) -
      3. * nz * (sqr(nx) + sqr(ny)) + cube(nx) + cube(ny)) {
    f = 0.; // solve case 3
  } else if (nx + ny <= nz) {
    f = nz + 0.5 * (nx + ny); // solve case 4
  } else {
    f = 0.; // solve case 5
  }

  return f - 0.5 * (nx + ny + nz);
}

  /*
  Scal a;
  Vect n(nx, ny, nz);

  Scal m1, m2, m3;
  m1 = nx;
  m2 = ny;
  m3 = nz;
  Scal m12 = m1 + m2;
  Scal pr = std::max(6.*m1*m2*m3, 1e-50);
  Scal V1 = m1*m1*m1/pr;
  Scal V2 = V1 + (m2 - m1)/(2.*m3), V3;
  Scal mm;
  if (m3 < m12) {
    mm = m3;
    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;
  }
  else {
    mm = m12;
    V3 = mm/(2.*m3);
  }

  Scal c = u;
  if (c < V1)
    a = pow (pr*c, 1./3.);
  else if (c < V2)
    a = (m1 + std::sqrt(m1*m1 + 8.*m2*m3*(c - V1)))/2.;
  else if (c < V3) {
    Scal p = 2.*m1*m2;
    Scal q = 3.*m1*m2*(m12 - 2.*m3*c)/2.;
    Scal p12 = std::sqrt(p);
    Scal teta = std::acos(q/(p*p12))/3.;
    Scal cs = std::cos(teta);
    a = p12*(std::sqrt(3.*(1. - cs*cs)) - cs) + m12;
  }
  else if (m12 < m3)
    a = m3*c + mm/2.;
  else {
    Scal p = m1*(m2 + m3) + m2*m3 - 1./4.;
    Scal q = 3.*m1*m2*m3*(1./2. - c)/2.;
    Scal p12 = std::sqrt(p);
    Scal teta = std::acos(q/(p*p12))/3.;
    Scal cs = std::cos(teta);
    a = p12*(std::sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }

  return a - (nx + ny + nz) * 0.5;
}
*/

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
  Scal nz = std::abs(n[2]);
  if (ny < nx) {
    std::swap(nx, ny);
  }
  if (nz < ny) {
    std::swap(ny, nz);
  }
  if (ny < nx) {
    std::swap(nx, ny);
  }
  
  Clip(u);

  if (u < 0.5) {
    return GetLineA0(nx, ny, nz, u);
  } else {
    return -GetLineA0(nx, ny, nz, 1. - u);
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

// GetLineU() helper
// assuming -0.5 * n.sum() < a < 0, 0 < nx < ny < nz
template <class Scal>
inline Scal GetLineU0(Scal nx, Scal ny, Scal nz, Scal a) {
  Scal f = 0.5 * (nx + ny + nz) + a;

  if (f <= 0) {
    return 0.;
  }

  if (nx > f) {
    return cube(f) / (6. * nx * ny * nz);
  } else if (ny > f) {
    return (3. * sqr(f) - 3. * f * nx + sqr(nx)) / (6. * ny * nz);
  } else if (nz > f && nx + ny > f) {
    nx = std::max(1e-50, nx);
    return (3. * sqr(f) - 3. * f * nx + sqr(nx) -
        std::min(1., (f - ny) / nx) * sqr(f - ny)) / (6. * ny * nz);
  } else if (nz >= f && nx + ny <= f) {
    return (2. * f - nx - ny) / (2. * nz);
  } else {
    return (cube(f) - cube(f - nx) - cube(f - ny) - cube(f - nz)) / 
      (6. * nx * ny * nz);
  }
}

// Sort to have a <= b <= c
template <class T>
inline void Sort(T& a, T& b, T& c) {
  if (b < a) {
    std::swap(a, b);
  }
  if (c < b) {
    std::swap(b, c);
  }
  if (b < a) {
    std::swap(a, b);
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
  Scal nz = std::abs(n[2]);

  Sort(nx, ny, nz);
  
  Clip(a, -0.5 * (nx + ny + nz), 0.5 * (nx + ny + nz));

  if (a < 0.) {
    return GetLineU0(nx, ny, nz, a);
  } else {
    return 1. - GetLineU0(nx, ny, nz, -a);
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

// Line ends by line constant
// n: normal
// a: line constant
// h: cell size
// Returns:
// two ends of segment inside cell (0,0 if no intersection)
template <class Scal>
inline std::array<GVect<Scal, 3>, 2> GetLineEnds(
    const GVect<Scal, 3>& n, Scal a, const GVect<Scal, 3>& h) {

  using Vect = GVect<Scal, 3>;
  // equation x.dot(n) = a;
  // (cell center is 0)
  Vect hh = h * 0.5;

  // intersection with -hh
  Vect xl((a + hh[1] * n[1]) / n[0], (a + hh[0] * n[0]) / n[1], 0); 
  // intersection with +hh
  Vect xr((a - hh[1] * n[1]) / n[0], (a - hh[0] * n[0]) / n[1], 0); 

  std::array<GVect<Scal, 3>, 2> e{Vect(0), Vect(0)}; // default to center
  size_t i = 0;

  if (-hh[0] <= xl[0] && xl[0] <= hh[0]) {
    e[i++] = Vect(xl[0], -hh[1], 0);
  } 
  if (-hh[0] <= xr[0] && xr[0] <= hh[0]) {
    e[i++] = Vect(xr[0], hh[1], 0);
  } 
  if (i < 2 && -hh[1] <= xl[1] && xl[1] <= hh[1]) {
    e[i++] = Vect(-hh[0], xl[1], 0);
  } 
  if (i < 2 && -hh[1] <= xr[1] && xr[1] <= hh[1]) {
    e[i++] = Vect(hh[0], xr[1], 0);
  } 
  if (i == 1) { // if only one point found, set second to the same
    e[i++] = e[0];
  } // if no points found, return default (cell center)
  return e;
}

// Line center by line constant
// n: normal
// a: line constant
// h: cell size
// Returns:
// mean point of intersection line (0 if no intersection)
template <class Scal>
inline GVect<Scal, 3> GetLineC(const GVect<Scal, 3>& n, Scal a,
                               const GVect<Scal, 3>& h) {
  using Vect = GVect<Scal, 3>;
  std::array<Vect, 2> e = GetLineEnds(n, a, h);
  return (e[0] + e[1]) * 0.5;
}

// Closest point to line
// x: target point
// n: normal
// a: line constant
// h: cell size
template <class Scal>
inline GVect<Scal, 3> GetNearest(const GVect<Scal, 3> x,
                                 const GVect<Scal, 3>& n, Scal a,
                                 const GVect<Scal, 3>& h) {
  using Vect = GVect<Scal, 3>;
  std::array<Vect, 2> e = GetLineEnds(n, a, h);
  Vect p = x + n * (e[0] - x).dot(n) / n.sqrnorm(); // projection to line
  if ((p - e[0]).dot(p - e[1]) < 0.) { // projection between ends
    return p;
  } else if ((x - e[0]).sqrnorm() < (x - e[1]).sqrnorm()) {
    return e[0];
  } 
  return e[1];
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
  MapFace<std::shared_ptr<CondFace>> mfc_;
  MapFace<std::shared_ptr<CondFace>> mfvz_; // zero-derivative bc for Vect

  FieldCell<Scal> fc_a_; // alpha (plane constant)
  FieldCell<Vect> fc_n_; // n (normal to plane)
  FieldCell<Scal> fc_us_; // smooth field
  FieldFace<Scal> ff_fu_; // volume flux
  FieldCell<Scal> fck_; // curvature
  FieldCell<Scal> fckp_; // curvature from particles
  FieldFace<Scal> ffu_; // field on faces
  FieldCell<Vect> fcg_; // gradient
  FieldFace<Vect> ffg_; // gradient
  size_t count_ = 0; // number of MakeIter() calls, used for splitting
  static constexpr size_t kNp = 5;
  static constexpr size_t kNpp = 1;
  FieldCell<std::array<Vect, kNp>> fcp_; // cell list
  FieldCell<std::array<Vect, kNp>> fcpt_; // cell list tmp
  FieldCell<std::array<Scal, kNp>> fcpw_; // cell list weight
  FieldCell<size_t> fcps_; // cell list size

 public:
  struct Par {
    bool curvgrad = false; // compute curvature using gradient
    bool part = false; // particles
    Scal part_relax = 1.; 
    Scal part_h0 = 1.; // dist init
    Scal part_h = 1.;  // dist eq
    Scal part_kstr = 1.; // stretching
    Scal part_kattr = 1.; // attraction to reconstructed interface
    Scal part_kbend = 1.; // bending
    size_t part_maxiter = 100; // num iter
    size_t part_dump_fr = 100; // num frames dump
    size_t part_report_fr = 100; // num frames report
    std::unique_ptr<Dumper> dmp; // dumper for particles
  };
  std::shared_ptr<Par> par;
  Par* GetPar() { return par.get(); }
  void SeedParticles(const FieldCell<Scal>& uc) {
    fcps_.Reinit(m, 0);
    fcp_.Reinit(m);
    fcpt_.Reinit(m);
    fcpw_.Reinit(m);

    auto& bc = m.GetBlockCells();
    auto& bn = m.GetBlockNodes();
    using MIdx = typename M::MIdx;
    MIdx wb = bn.GetBegin();
    Vect xb = m.GetNode(bn.GetIdx(wb));
    Vect h = m.GetNode(bn.GetIdx(wb + MIdx(1))) - xb;
    Scal hm = h.norminf();

    for (auto f : m.Faces()) {
      IdxCell cm = m.GetNeighbourCell(f, 0);
      IdxCell cp = m.GetNeighbourCell(f, 1);
      if (std::abs(uc[cm] - uc[cp]) > 1e-3 ) {
        //IdxCell c = (fc_n_[cp].norm() > fc_n_[cm].norm() ? cp : cm);
        IdxCell c = 
            (std::abs(uc[cp] - 0.5) < std::abs(uc[cm] - 0.5) ? cp : cm);
        Vect n = fc_n_[c];
        if (n.norm() > 1e-3) {
          n /= n.norm();
          //Vect x = m.GetCenter(c);
          Vect x = m.GetCenter(c) + GetLineC(fc_n_[c], fc_a_[c], h);
          Vect t = Vect(-n[1], n[0], 0.);
          if (fcps_[c] == 0) {
            for (int i = 0; i < kNp; ++i) {
              fcp_[c][fcps_[c]++] = 
                  x + t * (i - (kNp - 1) * 0.5) * hm * par->part_h0;
            }
          }
        }
      }
    }
  }
  void DumpParticles(size_t it) {
    // dump particles
    auto fr = par->part_dump_fr;
    size_t d = std::max<size_t>(1, par->part_maxiter / fr);
    if (fr > 1 && it % d == 0 || it + 1 == par->part_maxiter) {
      std::string st = "." + std::to_string(par->dmp->GetN());
      std::string sit = fr > 1 ? "_" + std::to_string(it) : "";
      std::string s = "partit" + st + sit + ".csv";
      std::cout 
          << "dump" 
          << " t=" << this->GetTime() + this->GetTimeStep()
          << " to " << s << std::endl;
      std::ofstream o;
      o.open(s);
      o << "x,y,z,c\n";

      for (auto c : m.Cells()) {
        for (size_t i = 0; i < fcps_[c]; ++i) {
          Vect x = fcp_[c][i];
          o << x[0] << "," << x[1] << "," << x[2] 
              << "," << (c.GetRaw() * 1234567 % 16) << "\n";
        }
      }
    }
  }
  Vof(M& m, const FieldCell<Scal>& fcu,
      const MapFace<std::shared_ptr<CondFace>>& mfc,
      const FieldFace<Scal>* ffv, const FieldCell<Scal>* fcs,
      double t, double dt, std::shared_ptr<Par> par)
      : AdvectionSolver<M>(t, dt, m, ffv, fcs)
      , mfc_(mfc), par(par)
      , fc_a_(m, 0), fc_n_(m, Vect(0)), fc_us_(m, 0), ff_fu_(m, 0) 
      , fck_(m, 0), fckp_(m, 0)
  {
    fc_u_.time_curr = fcu;
    for (auto it : mfc_) {
      IdxFace f = it.GetIdx();
      mfvz_[f] = std::make_shared<
          CondFaceGradFixed<Vect>>(Vect(0), it.GetValue()->GetNci());
    }
    Reconst(fc_u_.time_curr);
    if (par->part) {
      SeedParticles(fc_u_.time_curr);
      if (par->dmp->Try(this->GetTime() + this->GetTimeStep(), 
                        this->GetTimeStep())) {
        DumpParticles(par->part_maxiter - 1);
      }
    }

  }
  void StartStep() override {
    this->ClearIter();
    fc_u_.time_prev = fc_u_.time_curr;
    fc_u_.iter_curr = fc_u_.time_prev;
  }
  Scal Maxmod(Scal a, Scal b) {
    return std::abs(b) < std::abs(a) ? a : b;
  }
  // Normal with gradients
  void CalcNormal(const FieldCell<Scal>& uc) {
    auto uf = Interpolate(uc, mfc_, m);
    auto gc = Gradient(uf, m);
    for (auto c : m.AllCells()) {
      Vect g = gc[c];
      fc_n_[c] = g;
    }
  }
  // Normal with heigh function
  // XXX: no curvature
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
        if (d == Dir::i || std::abs(k) < std::abs(tk)) {
          tn = n;
          tk = k;
        } 
      }
      fc_n_[c] = tn;
    }
  }
  // Normal with heighfunction evaluated at only two points
  void CalcNormalHeightLite(const FieldCell<Scal>& uc) {
    using MIdx = typename M::MIdx;
    using Dir = typename M::Dir;
    auto& bc = m.GetBlockCells();
    for (auto c : m.SuCells()) {
      Vect tn; // bes[t] normal
      Scal tl; // bes[t] s[l]ope with minimal abs
      Scal tk; // bes[t] curvature[k]
      // direction of line tangent
      for (Dir d : {Dir::i, Dir::j}) {
        // direction of line normal ([d]irection [p]erpendicular)
        Dir dp(1 - size_t(d)); 

        MIdx w = bc.GetMIdx(c);

        // offset in dp
        MIdx op = MIdx(dp);

        // index shifted in d
        MIdx wm = bc.GetMIdx(c) - MIdx(d);
        MIdx wp = bc.GetMIdx(c) + MIdx(d);

        // height function 
        const Scal h = 
            uc[bc.GetIdx(w - op)] + 
            uc[bc.GetIdx(w)] + 
            uc[bc.GetIdx(w + op)];
        const Scal hm = 
            uc[bc.GetIdx(wm - op)] + 
            uc[bc.GetIdx(wm)] + 
            uc[bc.GetIdx(wm + op)];
        const Scal hp = 
            uc[bc.GetIdx(wp - op)] + 
            uc[bc.GetIdx(wp)] + 
            uc[bc.GetIdx(wp + op)];

        const Scal dx = m.GetCenter(c).dist(m.GetCenter(bc.GetIdx(wm)));

        
        // slope
        Scal lc = (hp - hm) * 0.5; // centered
        Scal l = lc; // XXX: choose centered approx
        // sign in dp
        Scal sg = uc[bc.GetIdx(w + MIdx(dp))] - uc[bc.GetIdx(w - MIdx(dp))];
        // curvature
        Scal k = -(hp - 2. * h + hm) / std::pow(1. + l * l, 3. / 2.) / dx;
        // normal
        Vect n;
        n[size_t(d)] = -l;
        n[size_t(dp)] = sg > 0. ? -1. : 1.;
        // select best with minimal abs
        if (d == Dir::i || std::abs(l) < std::abs(tl)) {
          tn = n;
          tl = l;
          tk = k;
        } 
      }
      fc_n_[c] = tn;
      Scal u = uc[c];
      const Scal th = 1e-6;
      fck_[c] = tk * (u > th && u < 1. - th ? 1. : 0.);
    }
  }
  void Reconst(const FieldCell<Scal>& uc) {
    //CalcNormal(uc);
    //CalcNormalHeight(uc);
    CalcNormalHeightLite(uc);
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
    //if (0) // XXX zero velocity
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

    if (par->curvgrad && sem("curv")) {
      ffu_ = Interpolate(fc_u_.iter_curr, mfc_, m); // [s]
      fcg_ = Gradient(ffu_, m); // [s]
      ffg_ = Interpolate(fcg_, mfvz_, m); // [i]

      fck_.Reinit(m); // curvature [i]
      for (auto c : m.Cells()) {
        Scal s = 0.;
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          auto& g = ffg_[f];
          // TODO: revise 1e-6
          auto n = g / (g.norm() + 1e-6);  // inner normal
          s += -n.dot(m.GetOutwardSurface(c, q));
        }
        fck_[c] = s / m.GetVolume(c);
      }
      m.Comm(&fck_);
    }

    if (par->part && sem("part")) {
      auto& uc = fc_u_.iter_curr;
      auto& bc = m.GetBlockCells();
      auto& bn = m.GetBlockNodes();
      using MIdx = typename M::MIdx;
      MIdx wb = bn.GetBegin();
      Vect xb = m.GetNode(bn.GetIdx(wb));
      Vect h = m.GetNode(bn.GetIdx(wb + MIdx(1))) - xb;
      Scal hm = h.norminf();

      SeedParticles(uc);

      const int sw = 1; // stencil width
      const int sn = sw * 2 + 1; // stencil size
      GBlock<IdxCell, dim> bo(MIdx(-sw, -sw, 0), MIdx(sn, sn, 1)); // offset

      bool dm = par->dmp->Try(this->GetTime() + this->GetTimeStep(), 
                              this->GetTimeStep());
      if (dm) {
        DumpParticles(0);
      }

      // advance particles
      // XXX: advance only if plan to dump
      if (dm) {
        for (size_t it = 0; it < par->part_maxiter; ++it) {
          // compute correction
          for (auto c : m.Cells()) {
            // clear force (needed for bending as applied from each corner)
            for (int i = 0; i < fcps_[c]; ++i) {
              fcpt_[c][i] = Vect(0); 
            }
            // traverse particles, append to force
            for (int i = 0; i < fcps_[c]; ++i) {
              const Vect& x = fcp_[c][i];
              Vect& t = fcpt_[c][i]; 
              // springs to adjacent particles
              for (int ii : {i - 1, i + 1}) {
                if (ii >= 0 && ii < fcps_[c]) {
                  const Vect& xx = fcp_[c][ii];
                  Vect dx = x - xx;
                  Scal d = par->part_h * hm - dx.norm();
                  t += dx / dx.norm() * d * par->part_kstr;
                }
              }

              // spring to nearest point at reconstructed interface
              auto w = bc.GetMIdx(c);
              bool fnd = false; // found at least one cell
              Vect xb;  // best point on the interface (nearest)
              for (auto wo : bo) {
                auto cc = bc.GetIdx(w + wo);
                if (!m.IsInner(cc)) {
                  continue;
                }
                Scal u = uc[cc];
                auto n = fc_n_[cc];
                if (u > 1e-3 && u < 1. - 1e-3) {
                  // nearest point to interface in cc
                  Vect xcc = m.GetCenter(cc);
                  Vect xn = xcc + GetNearest(x - xcc, n, fc_a_[cc], h);
                  if (!fnd || x.sqrdist(xn) < x.sqrdist(xb)) {
                    xb = xn;
                    fnd = true;
                  }
                }
              }
              if (fnd) {
                t += (xb - x) * par->part_kattr;
              }

              // bending
              if (i > 0 && i < fcps_[c] - 1) {
                int im = i - 1;
                int ip = i + 1;
                Vect xm = fcp_[c][im];
                Vect xp = fcp_[c][ip];
                Vect dm = xm - x;
                Vect dp = xp - x;
                // torque [length^2]
                Scal tq = par->part_kbend * 
                    (dm.norm() * dp.norm() + dm.dot(dp));
                // normal vectors [length]
                Vect nm(-dm[1], dm[0], 0.);
                Vect np(-dp[1], dp[0], 0.);
                // invert so they point inside angle
                nm *= (nm.dot(dp) > 0. ? 1. : -1.);
                np *= (np.dot(dm) > 0. ? 1. : -1.);
                // forces [length]
                Vect fm = nm * (-tq / dm.sqrnorm());
                Vect fp = np * (-tq / dp.sqrnorm());
                // apply
                fcpt_[c][im] += fm;
                fcpt_[c][ip] += fp;
                fcpt_[c][i] -= (fm + fp);
              }
            }
          }

          // compute correction norm
          Scal tmax = 0.;
          for (auto c : m.Cells()) {
            for (size_t i = 0; i < fcps_[c]; ++i) {
              fcp_[c][i] += fcpt_[c][i] * par->part_relax;
              tmax = std::max(tmax, fcpt_[c][i].norm());
            }
          }
          size_t dr = std::max<size_t>(1, 
              par->part_maxiter / par->part_report_fr);
          if (it % dr == 0 || it + 1 == par->part_maxiter) {
            std::cout << "it=" << it << " dxmax=" << tmax << std::endl;
          }

          // advance
          for (auto c : m.Cells()) {
            for (size_t i = 0; i < fcps_[c]; ++i) {
              fcp_[c][i] += fcpt_[c][i] * par->part_relax;
            }
          }

          if (dm) {
            DumpParticles(it + 1);
          }
        }

        // compute curvature
        fckp_.Reinit(m, 0.);
        for (auto c : m.Cells()) {
          // contains particles
          if (fcps_[c]) {
            int i = fcps_[c] / 2;
            int im = i - 1;
            int ip = i + 1;
            Vect x = fcp_[c][i];
            Vect xm = fcp_[c][im];
            Vect xp = fcp_[c][ip];
            Vect dm = xm - x;
            Vect dp = xp - x;
            Scal lm = dm.norm();
            Scal lp = dp.norm();
            Scal lmp = lm * lp;
            fckp_[c] = std::sqrt(2. * (lmp + dm.dot(dp))) / lmp;
          }
        }
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
  const FieldCell<Scal>& GetField(Layers l) const override {
    return fc_u_.Get(l);
  }
  const FieldCell<Scal>& GetAlpha() const {
    return fc_a_;
  }
  const FieldCell<Vect>& GetNormal() const {
    return fc_n_;
  }
  const FieldCell<Scal>& GetCurv() const override {
    return fck_;
  }
  // curvature from particles
  const FieldCell<Scal>& GetCurvP() const {
    return fckp_;
  }
  const FieldCell<std::array<Vect, kNp>>& GetPart() const {
    return fcp_;
  }
  const FieldCell<size_t>& GetPartS() const {
    return fcps_;
  }
  using P::GetField;
};

} // namespace solver
