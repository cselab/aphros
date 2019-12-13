#pragma once

#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "dump/vtk.h"
#include "func/init_u.h"
#include "func/primlist.h"
#include "hydro.h"
#include "parse/util.h"
#include "parse/vars.h"
#include "solver/advection.h"
#include "solver/approx.h"
#include "solver/cond.h"
#include "solver/fluid.h"

using namespace solver;
using namespace solver::fluid_condition;

template <class M>
FieldCell<typename M::Scal> GetBcField(MapCondFaceFluid& mf, const M& m) {
  FieldCell<typename M::Scal> fc(m, 0);
  for (auto it : mf) {
    IdxFace f = it.GetIdx();
    auto* b = it.GetValue().Get();
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

template <class M>
FieldCell<typename M::Vect> GetVort(
    const FieldCell<typename M::Vect>& fcv, const MapCondFace& mf, M& m) {
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
  } else if (vi == "grid_sin") {
    auto dx = var.Double["initvel_grid_sin_wavelength"];
    const Vect h = m.GetCellSize();
    for (auto c : m.AllCells()) {
      auto& v = fcv[c];
      auto x = m.GetCenter(c);
      v[0] = 0;
      v[1] = std::sin(x[0] / (dx * h[0]));
      v[2] = 0;
    }
  } else if (vi == "hill") {
    auto a = var.Double["hill_a"];
    auto b = var.Double["hill_b"];
    auto c = var.Double["hill_c"];
    Vect xc(var.Vect["hill_xc"]);
    for (auto i : m.AllCells()) {
      auto& v = fcv[i];
      auto xx = m.GetCenter(i) - xc;
      auto x = xx[0];
      auto y = xx[1];
      auto z = xx[2];
      Scal vrr; // vr / r
      Scal vz;
      if (x * x + y * y + z * z < 1.) { // inside sphere
        vrr = z;
        vz = -2 * pow(x, 2) - 2 * pow(y, 2) - pow(z, 2) + 1;
      } else { // outside
        vrr = z / pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 5.0 / 2.0);
        vz = (1.0 / 3.0) * (-pow(x, 2) - pow(y, 2) + 2 * pow(z, 2)) /
                 pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 5.0 / 2.0) -
             2.0 / 3.0;
      }

      v[0] = vrr * x;
      v[1] = vrr * y;
      v[2] = vz;

      // add swirl around z-axis
      v[0] += -2 * c * y / (x * x + y * y + b);
      v[1] += +2 * c * x / (x * x + y * y + b);

      v *= a;
    }
  } else if (vi == "hillorig") {
    auto a = var.Double["hill_a"];
    auto V = var.Double["hill_v"];
    // sphere r^2 + z^2 = a^2
    Vect xc(var.Vect["hill_xc"]);
    using std::pow;
    using std::sqrt;
    for (auto i : m.AllCells()) {
      auto& v = fcv[i];
      auto xx = m.GetCenter(i) - xc;
      auto x = xx[0];
      auto y = xx[1];
      auto z = xx[2];
      Scal vrr; // vr / r
      Scal vz;
      if (x * x + y * y + z * z < a * a) {
        vrr = (3.0 / 2.0) * V * z / pow(a, 2);
        vz = (1.0 / 2.0) * V *
             (5 * pow(a, 2) - 6 * pow(x, 2) - 6 * pow(y, 2) - 3 * pow(z, 2)) /
             pow(a, 2);
      } else {
        vrr = (3.0 / 2.0) * V * pow(a, 3) * z /
              pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 5.0 / 2.0);
        vz = (1.0 / 2.0) * V * pow(a, 3) *
             (-pow(x, 2) - pow(y, 2) + 2 * pow(z, 2)) /
             pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 5.0 / 2.0);
      }
      vz -= V;
      v[0] = vrr * x;
      v[1] = vrr * y;
      v[2] = vz;
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
    Vect n(var.Vect["ring_n"]); // normal
    n /= n.norm();
    Scal om = var.Double["ring_om"]; // vorticity in crosssection
    Scal r0 = var.Double["ring_r0"]; // inner radius
    Scal r1 = var.Double["ring_r1"]; // outer radius
    Scal qr = (r1 - r0) * 0.5; // radius
    Vect qc((r1 + r0) * 0.5, 0., 0.); // center
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
    Vect n(var.Vect["ring_n"]); // normal
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
          s += std::sin(iy * pi * y) * std::sin(iz * pi * z) /
               (iy * iz * (sqr(b) * sqr(iy) + sqr(iz)));
        }
      }
      fcv[i][0] = p * s * k;
    }
  } else if (vi == "uniform") {
    Vect v(var.Vect["vel"]);
    for (auto i : m.AllCells()) {
      fcv[i] = v;
    }
  } else if (vi == "box") {
    Vect v(var.Vect["vel"]);
    Vect a(var.Vect["vel_box_a"]);
    Vect b(var.Vect["vel_box_b"]);
    Rect<Vect> r(a, b);
    for (auto c : m.AllCells()) {
      if (r.IsInside(m.GetCenter(c))) {
        fcv[c] = v;
      }
    }
  } else if (vi == "solitonwang") {
    Scal xc(var.Double["soliton_xc"]);
    Scal yc(var.Double["soliton_yc"]);
    Scal yh(var.Double["soliton_yh"]);
    Scal xw(var.Double["soliton_xw"]);
    Scal e(var.Double["soliton_eps"]);
    using std::cos;
    using std::exp;
    using std::sin;
    using std::sqrt;
    Scal g(Vect(var.Vect["gravity"]).norm());
    Scal a = yh;
    Scal la = xw;
    Scal k = 2. * M_PI / la;
    for (auto c : m.AllCells()) {
      Vect xx = m.GetCenter(c);
      Scal x = xx[0] - xc;
      Scal y = xx[1];
      Scal kx = k * x;
      Scal ky = k * (y - yc);
      Scal h = yc + a / la *
                        (cos(kx) + 0.5 * e * cos(2 * kx) +
                         3. / 8 * sqr(e) * cos(3. * kx));
      Scal om = sqrt(g * k * (1 + sqr(e)));
      Scal u = om * a * exp(ky) * cos(kx);
      Scal v = om * a * exp(ky) * sin(kx);
      fcv[c] = Vect(u, v, 0.) * (y > h ? 0. : 1.);
    }
  } else if (vi == "solitoncos") {
    Scal xc(var.Double["soliton_xc"]);
    Scal yc(var.Double["soliton_yc"]);
    Scal yh(var.Double["soliton_yh"]);
    Scal xw(var.Double["soliton_xw"]);
    auto sinh = [](Scal t) { return std::sinh(t); };
    auto cosh = [](Scal t) { return std::cosh(t); };
    auto tanh = [](Scal t) { return std::tanh(t); };
    auto cos = [](Scal t) { return std::cos(t); };
    auto sin = [](Scal t) { return std::sin(t); };
    Scal g(Vect(var.Vect["gravity"]).norm());
    Scal p2 = 2 * M_PI;
    Scal L = xw;
    Scal D = yc;
    Scal C = std::sqrt(g * L * tanh(p2 * D / L) / p2);
    Scal T = L / C;
    Scal H = yh;
    for (auto c : m.AllCells()) {
      Vect xx = m.GetCenter(c);
      Scal x = xx[0] - xc;
      Scal y = xx[1];
      Scal Ux =
          H * p2 / T * cosh((p2 / L) * y) * cos(p2 * x / L) / sinh(p2 * D / L);
      Scal Uy =
          H * p2 / T * sinh((p2 / L) * y) * sin(p2 * x / L) / sinh(p2 * D / L);
      fcv[c] = Vect(Ux, Uy, 0.);
    }
  } else if (vi == "soliton") {
    Scal xc(var.Double["soliton_xc"]);
    Scal yc(var.Double["soliton_yc"]);
    Scal yh(var.Double["soliton_yh"]);

    using std::cosh;
    using std::sinh;
    using std::tanh;
    auto sech = [](Scal t) { return 1. / std::cosh(t); };
    using std::cos;
    using std::sin;
    using std::sqrt;

    Scal g(Vect(var.Vect["gravity"]).norm());
    Scal D = yc;
    Scal H = yh;

    H /= D;
    for (auto c : m.AllCells()) {
      Vect xx = m.GetCenter(c);
      Scal X = xx[0] - xc;
      Scal y = xx[1];

      X /= D;
      y /= D;
      Scal z = X * sqrt(3 * H / 4) * (1 - 5 * H / 8);
      Scal s = sech(z) * sech(z);
      Scal t = tanh(z);

      Scal u = H * (1 + H / 4 - 3 * H * y * y / 2) * s +
               H * H * (-1 + 9 * y * y / 4) * s * s;
      Scal v =
          sqrt(3) * pow(H, 3.0 / 2) * y * s * t *
          (1 - 3 * H / 8 - H * y * y / 2 + H * (-2 + 3 * y * y / 2) * s * s);
      Scal h = 1 + H * s - 4 * H * H * s * (1 - s) / 3;

      u *= sqrt(g * D);
      v *= sqrt(g * D);
      h *= D;
      y *= D;

      fcv[c] = Vect(u, v, 0.);
      if (y > h) {
        fcv[c] *= 0;
      }
    }
  } else if (vi == "wavelamb") {
    Scal a0(var.Double["wavelamb_a0"]);
    Scal xc(var.Double["wavelamb_xc"]);
    Scal h(var.Double["wavelamb_h"]);
    Scal k(var.Double["wavelamb_k"]);
    Scal g = -Vect(var.Vect["gravity"])[1];
    Vect air(var.Vect["wavelamb_air"]);
    Scal kvel(var.Double["wavelamb_kvel"]);

    using std::cos;
    using std::cosh;
    using std::pow;
    using std::sin;
    using std::sinh;
    using std::sqrt;
    using std::tanh;

    for (auto c : m.AllCells()) {
      auto xx = m.GetCenter(c);
      Scal x = xx[0] - xc;
      Scal y = xx[1] - h;
      Scal a = a0;

      Scal eps = a * k;
      Scal chi = 1.0 / tanh(h * k);
      Scal eta =
          (1.0 / 4.0) * a * chi * eps * (3 * pow(chi, 2) - 1) * cos(2 * k * x) +
          a * pow(eps, 2) *
              ((1.0 / 64.0) * (24 * pow(chi, 6) + 3 * pow(pow(chi, 2) - 1, 2)) *
                   cos(3 * k * x) +
               (1.0 / 8.0) * (-3 * pow(chi, 4) + 9 * pow(chi, 2) - 9) *
                   cos(k * x)) +
          a * cos(k * x);

      Scal omega = sqrt(
          g * k *
          (pow(eps, 2) * (pow(chi, 2) + (9.0 / 8.0) * pow(pow(chi, 2) - 1, 2)) +
           1) *
          tanh(h * k));

      Scal vx =
          (3.0 / 64.0) * a * pow(eps, 2) * g * k * (pow(chi, 2) - 1) *
              (pow(chi, 2) + 3) * (9 * pow(chi, 2) - 13) * cos(3 * k * x) *
              cosh(3 * k * (h + y)) / (omega * cosh(3 * h * k)) +
          a * g * k * cos(k * x) * cosh(k * (h + y)) / (omega * cosh(h * k)) +
          (3.0 / 4.0) * a * eps * g * k * pow(pow(chi, 2) - 1, 2) *
              cos(2 * k * x) * cosh(2 * k * (h + y)) /
              (chi * omega * cosh(2 * h * k));

      Scal vy =
          (3.0 / 64.0) * a * pow(eps, 2) * g * k * (pow(chi, 2) - 1) *
              (pow(chi, 2) + 3) * (9 * pow(chi, 2) - 13) * sin(3 * k * x) *
              sinh(3 * k * (h + y)) / (omega * cosh(3 * h * k)) +
          a * g * k * sin(k * x) * sinh(k * (h + y)) / (omega * cosh(h * k)) +
          (3.0 / 4.0) * a * eps * g * k * pow(pow(chi, 2) - 1, 2) *
              sin(2 * k * x) * sinh(2 * k * (h + y)) /
              (chi * omega * cosh(2 * h * k));

      fcv[c] = Vect(vx, vy, 0.) * kvel;
      if (y > eta) {
        fcv[c] = air;
      }
    }
  } else if (vi == "wavelamb_vort") {
    Scal a0(var.Double["wavelamb_a0"]);
    Scal xc(var.Double["wavelamb_xc"]);
    Scal h(var.Double["wavelamb_h"]);
    Scal k(var.Double["wavelamb_k"]);
    Scal d(var.Double["wavelamb_delta"]);
    Scal omk(var.Double["wavelamb_omk"]);
    Scal g = -Vect(var.Vect["gravity"])[1];
    Scal pi = M_PI;

    using std::cos;
    using std::cosh;
    using std::pow;
    using std::sin;
    using std::sinh;
    using std::sqrt;
    using std::tanh;

    for (auto c : m.AllCells()) {
      auto xx = m.GetCenter(c);
      Scal x = xx[0] - xc;
      Scal y = xx[1] - h;
      Scal a = a0;

      Scal eps = a * k;
      Scal chi = 1.0 / tanh(h * k);
      Scal eta =
          (1.0 / 4.0) * a * chi * eps * (3 * pow(chi, 2) - 1) * cos(2 * k * x) +
          a * pow(eps, 2) *
              ((1.0 / 64.0) * (24 * pow(chi, 6) + 3 * pow(pow(chi, 2) - 1, 2)) *
                   cos(3 * k * x) +
               (1.0 / 8.0) * (-3 * pow(chi, 4) + 9 * pow(chi, 2) - 9) *
                   cos(k * x)) +
          a * cos(k * x);

      Scal omega = sqrt(
          g * k *
          (pow(eps, 2) * (pow(chi, 2) + (9.0 / 8.0) * pow(pow(chi, 2) - 1, 2)) +
           1) *
          tanh(h * k));

      Scal vx =
          (3.0 / 64.0) * a * pow(eps, 2) * g * k * (pow(chi, 2) - 1) *
              (pow(chi, 2) + 3) * (9 * pow(chi, 2) - 13) * cos(3 * k * x) *
              cosh(3 * k * (h + y)) / (omega * cosh(3 * h * k)) +
          a * g * k * cos(k * x) * cosh(k * (h + y)) / (omega * cosh(h * k)) +
          (3.0 / 4.0) * a * eps * g * k * pow(pow(chi, 2) - 1, 2) *
              cos(2 * k * x) * cosh(2 * k * (h + y)) /
              (chi * omega * cosh(2 * h * k));

      Scal s = y - eta;

      Scal D = 1 / (sqrt(2 * pi * d * d)) * exp(-s * s / (2 * d * d));

      Scal omz = 2 * vx * D;

      fcv[c] = Vect(0., 0., omz) * omk;
    }
  } else if (vi == "solitonmccowan") {
    Scal xc(var.Double["soliton_xc"]);
    Scal yc(var.Double["soliton_yc"]);
    Scal yh(var.Double["soliton_yh"]);
    Scal N(var.Double["soliton_n"]);
    Scal MM(var.Double["soliton_m"]);

    using std::cosh;
    using std::sinh;
    using std::tanh;
    auto sech = [](Scal t) { return 1. / std::cosh(t); };
    using std::cos;
    using std::sin;
    using std::sqrt;

    Scal g(Vect(var.Vect["gravity"]).norm());
    Scal D = yc;
    Scal H = yh;

    H /= D;
    for (auto c : m.AllCells()) {
      Vect xx = m.GetCenter(c);
      Scal x = xx[0] - xc;
      Scal y = xx[1];

      x /= D;
      y /= D;
      Scal z = x * sqrt(3 * H / 4) * (1 - 5 * H / 8);
      Scal s = sech(z) * sech(z);
      Scal mx = MM * x;
      Scal my = MM * y;

      Scal C = 1. + 0.5 * H - 3. / 20. * sqr(H);
      Scal u = C * N * (1. + cos(my) * cosh(mx)) / sqr(cos(my) + cosh(mx));
      Scal v = C * N * (sin(my) * sinh(mx)) / sqr(cos(my) + cosh(mx));
      Scal h = 1 + H * s - 4 * H * H * s * (1 - s) / 3;

      u *= sqrt(g * D);
      v *= sqrt(g * D);
      h *= D;
      y *= D;

      fcv[c] = Vect(u, v, 0.);
      if (y > h) {
        fcv[c] *= 0;
      }
    }
  } else if (vi == "list") {
    auto pp = UPrimList<Scal>::ParseVel(
        var.String["vellist_path"], m.IsRoot(), var.Int["dim"]);
    fcv.Reinit(m, Vect(0));
    for (auto c : m.AllCells()) {
      Vect x = m.GetCenter(c);
      for (auto& p : pp) {
        fcv[c] += p.vel(x);
      }
    }
  } else if (vi == "zero") {
    // nop
  } else {
    throw std::runtime_error("Init(): unknown vel_init=" + vi);
  }
}

// argstr: argument string
// f: target face
// nc: target neighbour cell id
template <class M>
UniquePtr<CondFaceFluid> ParseFluidFaceCond(
    std::string argstr, IdxFace /*f*/, size_t nc, const M& /*m*/) {
  using namespace fluid_condition;
  using Vect = typename M::Vect;
  std::stringstream arg(argstr);

  std::string name;
  arg >> name;

  if (name == "wall") {
    // wall <velocity>
    // No-slip wall.
    // zero derivative for pressure, fixed for velocity,
    // fill-conditions for volume fraction.
    Vect vel;
    arg >> vel;
    return UniquePtr<NoSlipWallFixed<M>>(vel, nc);
  } else if (name == "inlet") {
    // inlet <velocity>
    // Fixed velocity inlet.
    Vect vel;
    arg >> vel;
    return UniquePtr<InletFixed<M>>(vel, nc);
  } else if (name == "inletflux") {
    // inletflux <velocity> <id>
    // Fixed flux inlet. Flux defined by given velocity is redistributed
    // over all faces with same id.
    Vect vel;
    int id;
    arg >> vel >> id;
    return UniquePtr<InletFlux<M>>(vel, id, nc);
  } else if (name == "outlet") {
    // Outlet. Velocity is extrapolated from neighbour cells and corrected
    // to yield zero total flux over outlet and inlet faces.
    return UniquePtr<OutletAuto<M>>(nc);
  } else if (name == "slipwall") {
    // Free-slip wall:
    // zero derivative for both pressure, velocity,
    // fill-conditions for volume fraction.
    // TODO: revise, should be non-penetration for velocity
    return UniquePtr<SlipWall<M>>(nc);
  } else if (name == "symm") {
    // Zero derivative for pressure, velocity and volume fraction
    // TODO: revise, should be non-penetration for velocity
    return UniquePtr<Symm<M>>(nc);
  }
  return nullptr;
}

template <class Scal>
std::ostream& operator<<(
    std::ostream& out, const CondFaceAdvection<Scal>& fca) {
  using Halo = typename CondFaceAdvection<Scal>::Halo;
  out << "nci=" << fca.nci << " clear0=" << fca.clear0
      << " clear1=" << fca.clear1
      << " halo=" << (fca.halo == Halo::fill ? "fill" : "reflect")
      << " fill_vf=" << fca.fill_vf << " fill_cl=" << fca.fill_cl;
  return out;
}

// Sets fields in cfa from arguments.
// Returns true if condition name is recognized.
// str: argument string
template <class Scal>
bool ParseAdvectionFaceCond(std::string str, CondFaceAdvection<Scal>& cfa) {
  using Halo = typename CondFaceAdvection<Scal>::Halo;
  std::stringstream arg(str);
  std::string name;
  arg >> name;

  if (name == "clear0") {
    Scal a;
    arg >> a;
    cfa.clear0 = a;
    return true;
  } else if (name == "clear1") {
    Scal a;
    arg >> a;
    cfa.clear1 = a;
    return true;
  } else if (name == "fill_vf") {
    Scal a;
    arg >> a;
    cfa.fill_vf = a;
    return true;
  } else if (name == "fill_cl") {
    Scal a;
    arg >> a;
    cfa.fill_cl = a;
    return true;
  } else if (name == "halo") {
    std::string halo;
    arg >> halo;
    if (halo == "fill") cfa.halo = Halo::fill;
    if (halo == "reflect") cfa.halo = Halo::reflect;
    return true;
  }
  return false;
}

template <class M>
void GetFluidFaceCond(
    const Vars& var, const M& m, MapCondFaceFluid& mff,
    MapCondFaceAdvection<typename M::Scal>& mfa) {
  using Dir = typename M::Dir;
  using MIdx = typename M::MIdx;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Halo = typename CondFaceAdvection<Scal>::Halo;
  size_t edim = var.Int["dim"];
  auto& fi = m.GetIndexFaces();
  MIdx gs = m.GetGlobalSize();

  // Returns true if face belongs to domain boundary
  auto gxm = [&fi](IdxFace i) -> bool {
    return fi.GetDir(i) == Dir::i && fi.GetMIdx(i)[0] == 0;
  };
  auto gxp = [&fi, gs](IdxFace i) -> bool {
    return fi.GetDir(i) == Dir::i && fi.GetMIdx(i)[0] == gs[0];
  };
  auto gym = [&fi](IdxFace i) -> bool {
    return fi.GetDir(i) == Dir::j && fi.GetMIdx(i)[1] == 0;
  };
  auto gyp = [&fi, gs](IdxFace i) -> bool {
    return fi.GetDir(i) == Dir::j && fi.GetMIdx(i)[1] == gs[1];
  };
  auto gzm = [&fi, edim](IdxFace i) -> bool {
    return edim >= 3 && fi.GetDir(i) == Dir::k && fi.GetMIdx(i)[2] == 0;
  };
  auto gzp = [&fi, gs, edim](IdxFace i) -> bool {
    return edim >= 3 && fi.GetDir(i) == Dir::k && fi.GetMIdx(i)[2] == gs[2];
  };

  using namespace solver::fluid_condition;
  using namespace solver;
  // default
  Scal clear0 = var.Double["bcc_clear0"];
  Scal clear1 = var.Double["bcc_clear1"];
  Scal inletcl = var.Double["inletcl"];
  Scal fill_vf = var.Double["bcc_fill"];

  auto set_adv = [&](IdxFace f, const std::vector<std::string>& ss) {
    static constexpr Scal kClNone = -1; // TODO define kClNone once
    auto& cb = mff[f];
    auto& cfa = mfa[f];
    cfa.nci = cb->GetNci();
    if (cb.Get<Symm<M>>()) {
      cfa.halo = Halo::reflect;
    } else if (cb.Get<NoSlipWall<M>>() || cb.Get<SlipWall<M>>()) {
      cfa.halo = Halo::fill;
      cfa.clear0 = clear0;
      cfa.clear1 = clear1;
      cfa.fill_vf = fill_vf;
      cfa.fill_cl = kClNone;
    } else if (cb.Get<Inlet<M>>()) {
      cfa.halo = Halo::fill;
      cfa.clear0 = 0;
      cfa.clear1 = 1;
      cfa.fill_vf = 0;
      cfa.fill_cl = inletcl;
    } else if (cb.Get<Outlet<M>>()) {
      cfa.halo = Halo::reflect;
    } else {
      throw std::runtime_error(std::string(__func__) + ": Unknown fluid bc");
    }

    for (auto s : ss) {
      if (!ParseAdvectionFaceCond(s, cfa)) {
        throw std::runtime_error("No advection condition found in '" + s + "'");
      }
    }
  };

  auto set_bc = [&](IdxFace f, size_t nci, std::string str) {
    auto& cff = mff[f];
    std::vector<std::string> ss; // strings not recognized as fluid cond
    for (auto s : Split(str, ',')) {
      if (!cff.Get<CondFaceFluid>()) {
        cff = ParseFluidFaceCond(s, f, nci, m);
        if (!cff.Get<CondFaceFluid>()) {
          ss.push_back(s);
        }
      } else {
        ss.push_back(s);
      }
    }
    if (!cff.Get<CondFaceFluid>()) {
      throw std::runtime_error("No fluid condition found in '" + str + "'");
    }
    set_adv(f, ss);
  };

  // List of domain boundaries: <name, func, nci>
  std::vector<std::tuple<std::string, std::function<bool(IdxFace)>, size_t>>
      bb = {{"bc_xm", gxm, 1}, {"bc_xp", gxp, 0}, {"bc_ym", gym, 1},
            {"bc_yp", gyp, 0}, {"bc_zm", gzm, 1}, {"bc_zp", gzp, 0}};

  // Set face conditions on domain boundaries
  for (auto& b : bb) {
    if (auto str = var.String(std::get<0>(b))) {
      for (auto f : m.AllFaces()) {
        if (std::get<1>(b)(f)) {
          set_bc(f, std::get<2>(b), *str);
        }
      }
    }
  }

  // selection boxes
  // Parameters (N>=0):
  // string boxN -- bc description
  // vect boxN_a -- lower corner
  // vect boxN_b -- upper corner
  // Check at least first nmax indices and all contiguous
  {
    int n = 0;
    const int nmax = 100;
    while (true) {
      std::string k = "box" + std::to_string(n);
      if (auto str = var.String(k)) {
        Vect a(var.Vect[k + "_a"]);
        Vect b(var.Vect[k + "_b"]);
        Rect<Vect> r(a, b);
        for (auto it : mff) {
          IdxFace f = it.GetIdx();
          if (r.IsInside(m.GetCenter(f))) {
            auto& cff = it.GetValue();
            auto nci = cff->GetNci();
            cff.Set(nullptr);
            set_bc(f, nci, *str);
          }
        }
      } else if (n > nmax) {
        break;
      }
      ++n;
    }
  }
  // inlet spheres
  // Parameters (N>=0):
  // vect inlet_sphN_c: center
  // vect inlet_sphN_r: radius
  // vect inlet_sphN_vel: maximum velocity
  // double inlet_sphN_vf: volume faction (0 or 1)
  // Check at least first nmax indices and all contiguous
  {
    int n = 0;
    const int nmax = 100;
    while (true) {
      std::string k = "inlet_sph" + std::to_string(n) + "_";
      if (auto p = var.Vect(k + "c")) {
        Vect xc(var.Vect[k + "c"]);
        Vect r(var.Vect[k + "r"]);
        Vect vel(var.Vect[k + "vel"]);
        Scal vf = var.Double[k + "vf"];

        // Returns intersection fraction of cell c and sphere
        auto V = [xc, r, &m](IdxCell c) {
          auto ls = [xc, r](const Vect& x) -> Scal {
            Vect xd = (x - xc) / r;
            return (1. - xd.sqrnorm()) * sqr(r.min());
          };
          return GetLevelSetVolume<Scal>(ls, m.GetCenter(c), m.GetCellSize());
        };

        for (auto it : mff) {
          IdxFace f = it.GetIdx();
          solver::CondFaceFluid* cb = it.GetValue().Get();
          auto nci = cb->GetNci();
          IdxCell c = m.GetNeighbourCell(f, nci);
          Scal inter = V(c);
          auto& cfa = mfa[f];
          if (inter > 0) {
            mff[f].Set<InletFixed<M>>(vel * inter, nci);
            cfa.halo = Halo::fill;
            cfa.clear0 = 0;
            cfa.clear1 = 1;
            cfa.fill_vf = vf;
            cfa.fill_cl = inletcl;
          }
        }
      } else if (n > nmax) {
        break;
      }
      ++n;
    }
  }
  /*
  // selection faces
  // Parameters (N>=0):
  // string faceN -- bc description
  // vect faceN_a -- lower corner
  // vect faceN_b -- upper corner
  // double faceN_vf -- inlet volume fraction
  // Check at least first nmax indices and all contiguous
  {
    int n = -1;
    const int nmax = 100;
    while (true) {
      ++n;
      std::string k = "face" + std::to_string(n);
      if (auto p = var.String(k)) {
        // set boundary conditions on faces of box (a,b)
        // normal to d with outer normal towards (b-a)[d]
        Vect a(var.Vect[k + "_a"]);
        Vect b(var.Vect[k + "_b"]);
        int d(var.Int[k + "_dir"]); // direction: 0:x, 1:y, 2:z
        Scal vf = var.Double[k + "_vf"];
        Rect<Vect> r(a, b);
        Vect h = m.GetCellSize();
        auto& cb = m.GetInBlockCells();
        auto& fi = m.GetIndexFaces();
        // indices of [a,b), [begin,end)
        Vect xd(0);
        xd[d] = 1.;
        // round to faces
        a = Vect(MIdx((a + h * 0.5) / h)) * h;
        b = Vect(MIdx((b + h * 0.5) / h)) * h;
        // indices
        MIdx wa(a / h + xd * 0.5);
        MIdx wb(b / h + xd * 0.5);
        wb[d] = wa[d] + 1; // size 1 in direction d
        // direction
        MIdx wd(0);
        wd[d] = 1;
        // direction of neighbour cell
        int nci = ((b - a)[d] > 0. ? 0 : 1);
        // box of valid indices
        MIdx w0 = cb.GetBegin();
        MIdx w1 = cb.GetEnd() + wd;
        // clip (a,b) to valid indices
        wa = wa.clip(w0, w1);
        wb = wb.clip(w0, w1);
        // size of local block
        MIdx ws = wb - wa;
        if (ws.prod() == 0) {
          continue;
        }
        typename M::BlockCells bb(wa, ws);
        for (auto w : bb) {
          IdxFace f = fi.GetIdx(w, Dir(d));
          mff[f] = ParseFluidFaceCond(*p, f, nci, m);
        }
        (void) vf;
      } else if (n > nmax) {
        break;
      }
    }
  }
  */
  // selection spheres
  // Parameters (N>=0):
  // string sphN -- bc description ("inlet")
  // vect sphN_c -- center
  // double sphN_r0 -- radius inner
  // double sphN_r1 -- radius outer
  // double sphN_vf -- inlet volume fraction
  // vect sphN_vel -- velocity
  // Check at least first nmax indices and all contiguous
  /*
  {
    int n = 0;
    const int nmax = 100;
    while (true) {
      std::string k = "sph" + std::to_string(n);
      if (auto p = var.String(k)) {
        if (*p == "inlet") {
          Vect xc(var.Vect[k + "_c"]);
          Scal r0 = var.Double[k + "_r0"];
          Scal r1 = var.Double[k + "_r1"];
          Scal vf = var.Double[k + "_vf"];
          Vect vel(var.Vect[k + "_vel"]);
          for (auto i : m.AllFaces()) {
            Vect x = m.GetCenter(i);
            if (xc.dist(x) < r1) {
              Scal r = xc.dist(x);
              Vect v = vel * std::min(1., (r1 - r) / (r1 - r0));
              std::string a = "inlet";
              a += " " + std::to_string(v[0]);
              a += " " + std::to_string(v[1]);
              a += " " + std::to_string(v[2]);
              if (set_bc(i, a)) {
                auto& b = mff[i];
                mfa[i] = UniquePtr<CondFaceValFixed<Scal>>(vf, b->GetNci());
              }
            }
          }
        } else {
          throw std::runtime_error("unknown selection sphere cond: " + *p);
        }
      } else if (n > nmax) {
        break;
      }
      ++n;
    }
  }
  */
}

template <class M>
void GetFluidCellCond(
    const Vars& var, M& m,
    MapCell<std::shared_ptr<solver::CondCellFluid>>& mcvel,
    std::pair<typename M::Scal, int>& pdist) {
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;

  auto sem = m.GetSem(__func__);
  if (sem("reduce")) {
    if (auto* p = var.Double("pfixed")) {
      Vect x(var.Vect["pfixed_x"]);
      // Find cell nearest to pfixed_x
      IdxCell c = m.FindNearestCell(x);
      pdist.first = m.GetCenter(c).dist(x);
      pdist.second = m.GetId();
      m.Reduce(std::make_shared<typename M::OpMinloc>(&pdist));
    }
  }

  if (sem("set")) {
    // Fix pressure at one cell
    if (auto* p = var.Double("pfixed")) {
      if (pdist.second == m.GetId()) {
        Vect x(var.Vect["pfixed_x"]);
        IdxCell c = m.FindNearestCell(x);
        mcvel[c] =
            std::make_shared<solver::fluid_condition::GivenPressureFixed<M>>(
                *p);
        std::cout << "pfixed id=" << pdist.second << " dist=" << pdist.first
                  << std::endl;
      }
    }

    // exclude cells
    int n = -1;
    const int nmax = 100;
    while (true) {
      ++n;
      std::string k = "cellbox" + std::to_string(n);
      if (auto p = var.String(k)) {
        // set boundary conditions on faces of box (a,b)
        // normal to d with outer normal towards (b-a)[d]
        Vect a(var.Vect[k + "_a"]);
        Vect b(var.Vect[k + "_b"]);
        Rect<Vect> r(a, b);
        Vect h = m.GetCellSize();
        auto& cb = m.GetInBlockCells();
        auto& ci = m.GetIndexCells();
        // round to faces
        a = Vect(MIdx((a + h * 0.5) / h)) * h;
        b = Vect(MIdx((b + h * 0.5) / h)) * h;
        // indices
        MIdx wa(a / h);
        MIdx wb(b / h);
        // box of valid indices
        MIdx w0 = cb.GetBegin();
        MIdx w1 = cb.GetEnd();
        // clip (a,b) to valid indices
        wa = wa.clip(w0, w1);
        wb = wb.clip(w0, w1);
        // size of local block
        MIdx ws = wb - wa;
        if (ws.prod() == 0) {
          continue;
        }
        typename M::BlockCells bb(wa, ws);
        for (auto w : bb) {
          IdxCell c = ci.GetIdx(w);
          mcvel[c] = std::make_shared<
              solver::fluid_condition::GivenVelocityAndPressureFixed<M>>(
              Vect(0), 0.);
        }
      } else if (n > nmax) {
        break;
      }
    }
  }
}

// Returns face polygon.
// x0,x1: points
// f0,f1: values
template <class M>
std::vector<typename M::Vect> GetPoly(IdxFace f, const M& m) {
  using Vect = typename M::Vect;
  std::vector<Vect> xx;
  for (size_t e = 0; e < m.GetNumNeighbourNodes(f); ++e) {
    auto n = m.GetNeighbourNode(f, e);
    xx.push_back(m.GetNode(n));
  }
  return xx;
}

template <class M>
void DumpBcFaces(
    const MapCondFaceAdvection<typename M::Scal>& mfa,
    const MapCondFaceFluid& mff, std::string fn, M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  auto sem = m.GetSem("DumpBcFaces");
  struct {
    std::vector<std::vector<Vect>> vxx;
    std::vector<Scal> vcond;
    std::vector<Scal> vcondf;
    std::vector<Scal> vblock;
  } * ctx(sem);
  auto& vxx = ctx->vxx;
  auto& vcond = ctx->vcond;
  auto& vcondf = ctx->vcondf;
  auto& vblock = ctx->vblock;
  if (sem("local")) {
    for (auto& it : mfa) {
      IdxFace f = it.GetIdx();
      vxx.push_back(GetPoly(f, m));
      const CondFaceAdvection<Scal>& b = it.GetValue();
      int cond = 0;
      int h = 0;
      using Halo = typename CondFaceAdvection<Scal>::Halo;
      switch (b.halo) {
        case Halo::reflect:
          h = 1;
          break;
        case Halo::fill:
          h = 2;
          break;
      }
      auto append = [&cond](int a) {
        a = std::min(99, std::max(0, a));
        cond = cond * 100 + a;
      };
      append(h);
      append(b.fill_vf * 10);
      append(b.fill_cl);
      append(b.clear0 * 10);
      append(b.clear1 * 10);

      Scal condf = -1;
      if (auto bs = mff.find(f)) {
        auto& b = *bs;
        if (b.Get<NoSlipWall<M>>()) {
          condf = 1;
        } else if (b.Get<SlipWall<M>>()) {
          condf = 2;
        } else if (b.Get<Inlet<M>>()) {
          condf = 3;
        } else if (b.Get<InletFlux<M>>()) {
          condf = 4;
        } else if (b.Get<Outlet<M>>()) {
          condf = 5;
        } else if (b.Get<Symm<M>>()) {
          condf = 6;
        }
      }
      vcond.push_back(cond);
      vcondf.push_back(condf);
      vblock.push_back(m.GetId());
    }

    using TV = typename M::template OpCatVT<Vect>;
    m.Reduce(std::make_shared<TV>(&vxx));
    using TS = typename M::template OpCatT<Scal>;
    m.Reduce(std::make_shared<TS>(&vcond));
    m.Reduce(std::make_shared<TS>(&vcondf));
    m.Reduce(std::make_shared<TS>(&vblock));
  }
  if (sem("write")) {
    if (m.IsRoot()) {
      std::cout << "dump"
                << " to " << fn << std::endl;
      WriteVtkPoly<Vect>(
          fn, vxx, nullptr, {&vcond, &vblock, &vcondf},
          {"advection", "block", "fluid"}, "Boundary conditions", true, true,
          true);
    }
  }
}
