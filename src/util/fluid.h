#pragma once

#include "solver/fluid.h"
#include "parse/vars.h"

// Returns field with the type (index)
// of boundary conditions in an adjacent face:
//   0: empty
//   1: no-slip wall
//   2: free-slip wall
//   3: inlet
//   4: outlet
//   -1: unknown
// mf: boundary conditions
template <class M>
FieldCell<typename M::Scal> GetBcField(
    MapFace<std::shared_ptr<solver::CondFaceFluid>>& mf, const M& m) {
  FieldCell<typename M::Scal> fc(m, 0);
  for (auto it : mf) {
    IdxFace f = it.GetIdx();
    auto* b = it.GetValue().get();
    size_t nci = b->GetNci();
    IdxCell c = m.GetNeighbourCell(f, nci);
    if (dynamic_cast<solver::fluid_condition::NoSlipWall<M>*>(b)) {
      fc[c] = 1.;
    } else if (dynamic_cast<solver::fluid_condition::SlipWall<M>*>(b)) {
      fc[c] = 2.;
    } else if (dynamic_cast<solver::fluid_condition::Inlet<M>*>(b)) {
      fc[c] = 3.;
    } else if (dynamic_cast<solver::fluid_condition::Outlet<M>*>(b)) {
      fc[c] = 4.;
    } else {
      fc[c] = -1.;
    }
  }
  return fc;
}

// Computes vorticity of vector field.
// fcv: vector field [s]
// mf: boundary conditions for fcv
// Returns:
// fco: vorticity [i]
template <class M>
FieldCell<typename M::Vect> GetVort(const FieldCell<typename M::Vect>& fcv,
                       const MapFace<std::shared_ptr<solver::CondFace>>& mf,
                       M& m) {
  auto ffv = solver::Interpolate(fcv, mf, m);

  auto d0 = solver::Gradient(GetComponent(ffv, 0), m);
  auto d1 = solver::Gradient(GetComponent(ffv, 1), m);
  auto d2 = solver::Gradient(GetComponent(ffv, 2), m);

  FieldCell<typename M::Vect> r(m);
  for (auto c : m.Cells()) {
    r[c][0] = d2[c][1] - d1[c][2];
    r[c][1] = d0[c][2] - d2[c][0];
    r[c][2] = d1[c][0] - d0[c][1];
  }

  return r;
}


// Initializes velocity from parameters.
// fcv: vector field [s]
// mf: boundary conditions for fcv
// Returns:
// fco: vorticity [i]
template <class M>
void InitVel(FieldCell<typename M::Vect>& fcv, const Vars& var, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  const std::string vi = var.String["vel_init"];
  if (vi == "taylor-green") {
    for (auto i : m.AllCells()) {
      auto& v = fcv[i];
      auto x = m.GetCenter(i);
      if (var.Int["dim"] == 2) {
        x[2] = 0.;
      }
      v[0] = std::sin(x[0]) * std::cos(x[1]) * std::cos(x[2]);
      v[1] = -std::cos(x[0]) * std::sin(x[1]) * std::cos(x[2]);
      v[2] = 0.;
    }
  } else if (vi == "vortex") {
    Scal pi = M_PI;
    Vect xc1(var.Vect["vort_x1"]);
    Vect xc2(var.Vect["vort_x2"]);
    Scal g1 = var.Double["vort_g1"];
    Scal g2 = var.Double["vort_g2"];
    Scal vmax = var.Double["vort_vmax"];
    std::vector<Vect> xxc = {xc1, xc2};
    std::vector<Scal> gg = {g1, g2};
    for (size_t i = 0; i < xxc.size(); ++i) {
      Vect xc = xxc[i];
      Scal g = gg[i];
      for (auto c : m.AllCells()) {
        Vect x = m.GetCenter(c);
        Scal r = xc.dist(x);
        Scal r2 = r * r;
        Vect v(0);
        v[0] = -(x[1] - xc[1]) * g / r2 * 2. * pi;
        v[1] = (x[0] - xc[0]) * g / r2 * 2. * pi;
        if (v.norm() > vmax) {
          v *= vmax / v.norm();
        }
        fcv[c][0] += v[0];
        fcv[c][1] += v[1];
      }
    }
  } else if (vi == "vortexring") { // XXX: use with initvort=1
    Vect xc(var.Vect["ring_c"]); // center of ring
    Vect n(var.Vect["ring_n"]);  // normal
    n /= n.norm();
    Scal om = var.Double["ring_om"];  // vorticity in crosssection
    Scal r0 = var.Double["ring_r0"];  // inner radius
    Scal r1 = var.Double["ring_r1"];  // outer radius
    Scal qr = (r1 - r0) * 0.5;   // radius
    Vect qc((r1 + r0) * 0.5, 0., 0.);   // center
    const Scal eps = 1e-10;
    for (auto c : m.AllCells()) {
      Vect x = m.GetCenter(c) - xc;
      // along axis
      Scal xn = n.dot(x);
      // along plane
      Scal xt = (x - n * xn).norm();
      // unit along plane
      Vect t = (x - n * xn) / std::max(eps, xt);
      // unit along circle 
      Vect s = n.cross(t);
      Vect q(xt, xn, 0.);
      fcv[c] = ((q - qc).sqrnorm() <= sqr(qr) ? s * om : Vect(0));
    }
  } else if (vi == "vortexgauss") {
    // XXX: Bergdorf 2007, Direct numerical simulations of vortex rings
    Scal pi = M_PI;
    Scal g(var.Double["ring_gamma"]); // circulation
    Scal sig(var.Double["ring_sigma"]); // sigma
    Scal rad(var.Double["ring_r"]); // radius
    Scal nfr(var.Double["ring_noise_freq"]); // noise angular frequency
    Scal namp(var.Double["ring_noise_amp"]); // noise amp relative to r
    Scal nfr2(var.Double["ring_noise2_freq"]); // noise angular frequency
    Scal namp2(var.Double["ring_noise2_amp"]); // noise amp relative to r
    Vect xc(var.Vect["ring_c"]); // center
    Vect n(var.Vect["ring_n"]);  // normal
    n /= n.norm();
    const Scal eps = 1e-10;
    for (auto c : m.AllCells()) {
      Vect x = m.GetCenter(c) - xc;
      // along axis
      Scal xn = n.dot(x);
      // select direction 
      Vect v(0);
      v[n.abs().argmin()] = 1.;
      // unit vectors in plane
      Vect vx = n.cross(v);
      Vect vy = n.cross(vx);
      // angle
      Scal a = std::atan2(vx.dot(x), vy.dot(x));
      // noise
      xn += rad * namp * std::sin(a * nfr);
      // along plane
      Scal xt = (x - n * xn).norm();
      // noise2
      xt += rad * namp2 * std::sin(a * nfr2);
      // unit radial along plane
      Vect et = (x - n * xn) / std::max(eps, xt);
      // unit along circle 
      Vect es = n.cross(et);
      Scal s2 = sqr(xn) + sqr(xt - rad);
      Scal sig2 = sqr(sig);
      Scal om = g / (pi * sig2) * std::exp(-s2 / sig2);
      fcv[c] = es * om;
    }
  } else if (vi == "rot") {
    Vect xc(var.Vect["rot_c"]);
    Scal om = var.Double["rot_om"];
    for (auto c : m.AllCells()) {
      Vect x = m.GetCenter(c);
      Vect v(0);
      v[0] = -(x[1] - xc[1]) * om * 0.5;
      v[1] = (x[0] - xc[0]) * om * 0.5;
      fcv[c] = v;
    }
  } else if (vi == "grad") {
    Vect xc(var.Vect["grad_c"]);
    Vect vx(var.Vect["grad_vx"]); // gradient of vel[0]
    Vect vy(var.Vect["grad_vy"]); // gradient of vel[1]
    Vect vz(var.Vect["grad_vz"]); // gradient of vel[2]
    for (auto c : m.AllCells()) {
      Vect x = m.GetCenter(c) - xc;
      Vect v;
      v[0] = vx.dot(x);
      v[1] = vy.dot(x);
      v[2] = vz.dot(x);
      fcv[c] = v;
    }
  } else if (vi == "pois" || vi == "poisy") {
    // Poiseuille with walls in y
    auto gs = m.GetGlobalSize(); // global mesh size
    Scal ext = var.Double["extent"]; // TODO: revise
    Vect gh = Vect(gs) * ext / gs.max(); // global domain length
    Scal pv = var.Double["poisvel"]; // centerline velocity

    for (auto i : m.AllCells()) {
      Scal y = m.GetCenter(i)[1] / gh[1];
      fcv[i][0] = y * (1. - y) * 4. * pv;
    }
  } else if (vi == "poisyz") {
    // Poiseuille with walls in y and z
    // Spiga 1994: Symmetric solution for velocity in rectangular ducts
    auto gs = m.GetGlobalSize(); // global mesh size
    Scal ext = var.Double["extent"]; // TODO: revise
    Vect gh = Vect(gs) * ext / gs.max(); // global domain length
    Scal mu = var.Double["poismu"]; // viscosity
    Scal pg = var.Double["poisgrad"]; // pressure gradient
    int im = var.Int["poisiter"]; // depth to evaluate series
    bool wym = var.Int["poiswym"]; // wallym
    bool wyp = var.Int["poiswyp"]; // wallyp
    bool wzm = var.Int["poiswzm"]; // wallzm
    bool wzp = var.Int["poiswzp"]; // wallzp
    Scal pi = M_PI;

    Scal ly = gh[1];
    Scal lz = gh[2];

    if ((!wym && !wyp) || (!wzm && !wzp)) {
      throw std::runtime_error("poisyz: can't remove both walls");
    }

    Scal oy = 0.;
    if (!wym) {
      oy = 0.5;
      ly *= 2;
    } else if (!wyp) {
      oy = 0;
      ly *= 2;
    }

    Scal oz = 0.;
    if (!wzm) {
      oz = 0.5;
      lz *= 2;
    } else if (!wzp) {
      oz = 0;
      lz *= 2;
    }


    Scal p = sqr(ly) * pg / mu;
    Scal b = lz / ly;
    Scal k = 16. * sqr(b) / std::pow(pi, 4);

    // TODO: tests for gh[1] != gh[2]
    for (auto i : m.AllCells()) {
      Scal y = m.GetCenter(i)[1] / ly;
      y = y + oy;
      Scal z = m.GetCenter(i)[2] / lz;
      z = z + oz;
      Scal s = 0.;
      for (int iy = 1; iy < im * 2; iy += 2) {
        for (int iz = 1; iz < im * 2; iz += 2) {
          s += std::sin(iy * pi * y) *  std::sin(iz * pi * z) / 
              (iy * iz * (sqr(b) * sqr(iy) + sqr(iz)));
        }
      }
      fcv[i][0] = p * s * k;
    }
  } else if (vi == "uniform" ) {
    Vect v(var.Vect["vel"]);
    for (auto i : m.AllCells()) {
      fcv[i] = v;
    }
  } else if (vi == "zero" ) {
    // nop
  } else  {
    throw std::runtime_error("Init(): unknown vel_init=" + vi);
  }
}

