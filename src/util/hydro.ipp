// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <functional>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "convdiff.h"
#include "dump/dumper.h"
#include "dump/vtk.h"
#include "func/init_u.h"
#include "func/primlist.h"
#include "parse/util.h"
#include "solver/approx.h"
#include "solver/cond.h"
#include "solver/pois.h"
#include "solver/sphavg.h"
#include "solver/vof.h"
#include "solver/vof_eb.h"
#include "solver/vofm.h"

#include "hydro.h"

using namespace fluid_condition;

template <class M>
FieldCell<typename M::Scal> GetBcField(MapCondFaceFluid& mf, const M& m) {
  FieldCell<typename M::Scal> fc(m, 0);
  for (const auto& it : mf) {
    const IdxFace f = it.first;
    auto& cb = it.second;
    const size_t nci = cb->GetNci();
    const IdxCell c = m.GetCell(f, nci);
    if (cb.template Get<fluid_condition::NoSlipWall<M>>()) {
      fc[c] = 1.;
    } else if (cb.template Get<fluid_condition::SlipWall<M>>()) {
      fc[c] = 2.;
    } else if (cb.template Get<fluid_condition::Inlet<M>>()) {
      fc[c] = 3.;
    } else if (cb.template Get<fluid_condition::Outlet<M>>()) {
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
  auto ffv = Interpolate(fcv, mf, m);

  auto d0 = Gradient(GetComponent(ffv, 0), m);
  auto d1 = Gradient(GetComponent(ffv, 1), m);
  auto d2 = Gradient(GetComponent(ffv, 2), m);

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

// Initializes advection condition
// with default values given a fluid condition.
template <class M, class Scal = typename M::Scal>
void SetDefaultAdvectionFaceCond(
    CondFaceAdvection<Scal>& ca, const UniquePtr<CondFaceFluid>& cf,
    Scal clear0, Scal clear1, Scal inletcl, Scal fill_vf) {
  using Halo = typename CondFaceAdvection<Scal>::Halo;
  static constexpr Scal kClNone = -1; // TODO define kClNone once
  ca.nci = cf->GetNci();
  if (cf.Get<Symm<M>>()) {
    ca.halo = Halo::reflect;
  } else if (cf.Get<NoSlipWall<M>>() || cf.Get<SlipWall<M>>()) {
    ca.halo = Halo::fill;
    ca.clear0 = clear0;
    ca.clear1 = clear1;
    ca.fill_vf = fill_vf;
    ca.fill_cl = kClNone;
  } else if (cf.Get<Inlet<M>>()) {
    ca.halo = Halo::fill;
    ca.clear0 = 0;
    ca.clear1 = 1;
    ca.fill_vf = 0;
    ca.fill_cl = inletcl;
  } else if (cf.Get<Outlet<M>>()) {
    ca.halo = Halo::reflect;
  } else {
    throw std::runtime_error(std::string(__func__) + ": Unknown fluid bc");
  }
}

template <class M, class Scal = typename M::Scal>
void ParseFaceCond(
    IdxFace f, size_t nci, std::string str, const M& m, Scal clear0,
    Scal clear1, Scal inletcl, Scal fill_vf, UniquePtr<CondFaceFluid>& cf,
    CondFaceAdvection<typename M::Scal>& ca) {
  std::vector<std::string> ss; // strings not recognized as fluid cond
  for (auto s : Split(str, ',')) {
    auto cft = ParseFluidFaceCond(s, f, nci, m); // try to parse as fluid
    if (!cft.Get()) { // otherwise try as advection later
      ss.push_back(s);
    } else {
      cf = std::move(cft);
    }
  }
  if (!cf.Get<CondFaceFluid>()) {
    throw std::runtime_error("No fluid condition found in '" + str + "'");
  }
  SetDefaultAdvectionFaceCond<M>(ca, cf, clear0, clear1, inletcl, fill_vf);
  for (auto s : ss) {
    if (!ParseAdvectionFaceCond(s, ca)) {
      throw std::runtime_error("No advection condition found in '" + s + "'");
    }
  }
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

  using namespace fluid_condition;
  // default
  const Scal clear0 = var.Double["bcc_clear0"];
  const Scal clear1 = var.Double["bcc_clear1"];
  const Scal inletcl = var.Double["inletcl"];
  const Scal fill_vf = var.Double["bcc_fill"];

  auto set_bc = [&](IdxFace f, size_t nci, std::string str) {
    ParseFaceCond(
        f, nci, str, m, clear0, clear1, inletcl, fill_vf, mff[f], mfa[f]);
  };

  // List of domain boundaries: <name, func, nci>
  std::vector<std::tuple<std::string, std::function<bool(IdxFace)>, size_t>>
      bb = {{"bc_xm", gxm, 1}, {"bc_xp", gxp, 0}, {"bc_ym", gym, 1},
            {"bc_yp", gyp, 0}, {"bc_zm", gzm, 1}, {"bc_zp", gzp, 0}};

  // Set face conditions on domain boundaries
  for (auto& b : bb) {
    if (auto str = var.String.Find(std::get<0>(b))) {
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
      if (auto str = var.String.Find(k)) {
        Vect a(var.Vect[k + "_a"]);
        Vect b(var.Vect[k + "_b"]);
        Rect<Vect> r(a, b);
        for (auto& it : mff) {
          const IdxFace f = it.first;
          if (r.IsInside(m.GetCenter(f))) {
            auto& cff = it.second;
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
      if (auto p = var.Vect.Find(k + "c")) {
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

        for (auto& it : mff) {
          IdxFace f = it.first;
          CondFaceFluid* cb = it.second.Get();
          auto nci = cb->GetNci();
          IdxCell c = m.GetCell(f, nci);
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
      if (auto str = var.String.Find(k)) {
        // set boundary conditions on faces of box (a,b)
        // normal to d with outer normal towards (b-a)[d]
        Vect a(var.Vect[k + "_a"]);
        Vect b(var.Vect[k + "_b"]);
        int d(var.Int[k + "_dir"]); // direction: 0:x, 1:y, 2:z
        Rect<Vect> r(a, b);
        Vect h = m.GetCellSize();
        auto& cb = m.GetAllBlockCells();
        auto& fi = m.GetIndexFaces();
        Vect xd(0);
        xd[d] = 1.;
        // round to faces
        a = Vect(MIdx((a + h * 0.5) / h)) * h;
        b = Vect(MIdx((b + h * 0.5) / h)) * h;
        // indices of [a,b), [begin,end)
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
        if ((wb - wa).prod() == 0) {
          continue;
        }
        typename M::BlockCells bb(wa, wb - wa);
        for (auto w : bb) {
          IdxFace f = fi.GetIdx(w, Dir(d));
          set_bc(f, nci, *str);
        }
      } else if (n > nmax) {
        break;
      }
    }
  }
}

template <class M>
void GetFluidCellCond(
    const Vars& var, M& m, MapCell<std::shared_ptr<CondCellFluid>>& mcvel) {
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;

  auto sem = m.GetSem(__func__);

  struct {
    std::pair<typename M::Scal, int> pdist;
  } * ctx(sem);
  auto& pdist = ctx->pdist;

  if (sem("reduce")) {
    if (auto* p = var.Double.Find("pfixed")) {
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
    if (auto* p = var.Double.Find("pfixed")) {
      if (pdist.second == m.GetId()) {
        Vect x(var.Vect["pfixed_x"]);
        IdxCell c = m.FindNearestCell(x);
        mcvel[c] = std::make_shared<fluid_condition::GivenPressureFixed<M>>(*p);
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
      if (auto p = var.String.Find(k)) {
        // set boundary conditions on faces of box (a,b)
        // normal to d with outer normal towards (b-a)[d]
        Vect a(var.Vect[k + "_a"]);
        Vect b(var.Vect[k + "_b"]);
        Rect<Vect> r(a, b);
        Vect h = m.GetCellSize();
        auto& cb = m.GetAllBlockCells();
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
              fluid_condition::GivenVelocityAndPressureFixed<M>>(Vect(0), 0.);
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
  for (size_t e = 0; e < m.GetNumNodes(f); ++e) {
    auto n = m.GetNode(f, e);
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
      IdxFace f = it.first;
      const CondFaceAdvection<Scal>& b = it.second;
      if (!m.IsInner(m.GetCell(f, b.GetNci()))) {
        continue;
      }
      vxx.push_back(GetPoly(f, m));
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

// FIXME: only fills SuFaces, more halo faces AllFaces may be needed
template <class M, class Scal>
void AppendBodyCond(
    const FieldCell<bool>& fc, std::string str, const M& m, Scal clear0,
    Scal clear1, Scal inletcl, Scal fill_vf,
    MapCell<std::shared_ptr<CondCellFluid>>* mcf, MapCondFaceFluid& mff,
    MapCondFaceAdvection<Scal>& mfa) {
  using Vect = typename M::Vect;
  for (auto f : m.SuFaces()) {
    const IdxCell cm = m.GetCell(f, 0);
    const IdxCell cp = m.GetCell(f, 1);
    if (fc[cm] != fc[cp]) {
      size_t nci = (fc[cm] ? 1 : 0); // inner cell
      auto& condf = mff[f];
      auto& conda = mfa[f];
      ParseFaceCond(
          f, nci, str, m, clear0, clear1, inletcl, fill_vf, condf, conda);
      if (mcf) {
        (*mcf)[m.GetCell(f, 1 - nci)] =
            std::make_shared<fluid_condition::GivenVelocityAndPressureFixed<M>>(
                Vect(0), 0.);
      }
    }
  }
}

// Computes velocity fcvel from vorticity fcvort
template <class M, class Vect>
void InitVort(
    const FieldCell<Vect>& fcvort, FieldCell<Vect>& fcvel,
    const MapCondFaceFluid& mf_fluid, bool verb, M& m) {
  using Scal = typename M::Scal;
  auto sem = m.GetSem("initvort");
  struct {
    MapCondFace mfcw; // velocity cond
    MapCondFace mfcwd; // velocity cond in one diretion
    std::shared_ptr<PoisSolver<M>> ps;
    FieldCell<Scal> fcvorts; // one component of vorticity
    FieldCell<Vect> fcpot; // velocity potential
  } * ctx(sem);
  if (sem("initpois")) {
    ctx->fcpot.Reinit(m);
    ctx->mfcw = GetVelCond(m, mf_fluid);
  }
  for (size_t d = 0; d < M::dim; ++d) {
    const std::string dname = std::to_string(d);
    if (sem("init-" + dname)) {
      ctx->mfcwd = GetScalarCond(ctx->mfcw, d, m);
      ctx->ps = std::make_shared<PoisSolver<M>>(ctx->mfcwd, m);
      ctx->fcvorts = GetComponent(fcvort, d);
      for (auto c : m.Cells()) {
        ctx->fcvorts[c] *= -1;
      }
    }
    if (sem.Nested("solve-" + dname)) {
      ctx->ps->Solve(ctx->fcvorts);
    }
    if (sem("post-" + dname)) {
      SetComponent(ctx->fcpot, d, ctx->ps->GetField());
      if (m.IsRoot() && verb) {
        std::cout << "om" << ("xyz"[d]) << ":"
                  << " res=" << m.GetResidual() << " iter=" << m.GetIter()
                  << std::endl;
      }
    }
  }
  if (sem("vel")) {
    fcvel = GetVort(ctx->fcpot, ctx->mfcw, m);
    m.Comm(&fcvel);
  }
}

template <class M>
void DumpTraj(
    M& m, bool dm, const Vars& var, size_t frame, typename M::Scal t,
    const GRange<size_t>& layers,
    const Multi<const FieldCell<typename M::Scal>*>& fcvf,
    const Multi<const FieldCell<typename M::Scal>*>& fccl,
    const Multi<const FieldCell<typename M::MIdx>*>& fcim,
    const FieldCell<typename M::Scal>& fcp,
    const FieldCell<typename M::Vect>& fcvel,
    const FieldCell<typename M::Vect>& fcvelm, typename M::Scal dt) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using SA = Sphavg<M>; // spherical averages
  constexpr Scal kClNone = -1;

  auto sem = m.GetSem("dumptraj");
  struct {
    std::vector<std::string> names; // color reduce: variable name
    std::vector<Scal> colors; // color reduce: cl
    std::vector<std::vector<Scal>> values; // color reduce: vector
    std::unique_ptr<SA> sphavg; // spherical averages
    std::vector<typename SA::Sph> vsph;
  } * ctx(sem);
  auto& names = ctx->names;
  auto& colors = ctx->colors;
  auto& values = ctx->values;
  auto& sphavg = ctx->sphavg;
  auto& vsph = ctx->vsph;
  if (sem("color-calc")) {
    std::map<Scal, std::vector<Scal>> mp; // map color to vector
    Vect gh = m.GetGlobalLength(); // global domain length

    // add scalar name
    auto nma = [&](const std::string nm) { names.push_back(nm); };
    // add vector name
    auto nmav = [&](const std::string nm) {
      names.push_back(nm + "x");
      names.push_back(nm + "y");
      names.push_back(nm + "z");
    };

    // XXX: adhoc, the following order assumed in post: vs,r,x,y,z,...

    // list of vars // TODO: revise
    names.clear();
    nma("vf");
    nma("r");
    nmav("");
    nmav("v");
    nma("p");
    nma("xx");
    nma("xy");
    nma("xz");
    nma("yy");
    nma("yz");
    nma("zz");

    // traverse cells, append to mp
    for (auto c : m.Cells()) {
      for (auto l : layers) {
        auto cl = (*fccl[l])[c];
        if (cl != kClNone) {
          auto& v = mp[cl]; // vector for data
          auto x = m.GetCenter(c); // cell center
          x += Vect((*fcim[l])[c]) * gh; // translation by image vector
          auto w = (*fcvf[l])[c] * m.GetVolume(c); // volume

          size_t i = 0;
          // append scalar value
          auto add = [&v, &i](Scal a) {
            if (i >= v.size()) {
              v.resize(i + 1);
            }
            v[i] += a;
            ++i;
          };
          // append vector value
          auto addv = [&](Vect a) {
            add(a[0]);
            add(a[1]);
            add(a[2]);
          };

          // list of vars, XXX: keep consistent with names
          add(w); // vf,  XXX: adhoc, vf must be first, divided on dump
          add(0.); // r,  XXX: adhoc, r must be second, computed on dump
          addv(x * w); // x
          addv(fcvel[c] * w); // velocity
          add(fcp[c] * w); // pressure
          add(x[0] * x[0] * w); // xx
          add(x[0] * x[1] * w); // xy
          add(x[0] * x[2] * w); // xz
          add(x[1] * x[1] * w); // yy
          add(x[1] * x[2] * w); // yz
          add(x[2] * x[2] * w); // zz
        }
      }
    }
    // copy to vector
    colors.clear();
    values.clear();
    for (auto& it : mp) {
      colors.push_back(it.first); // color
      values.push_back(it.second); // vector
    }
    using TS = typename M::template OpCatT<Scal>;
    using TVS = typename M::template OpCatVT<Scal>;
    m.Reduce(std::make_shared<TS>(&colors));
    m.Reduce(std::make_shared<TVS>(&values));
  }
  if (sem("color-post")) {
    if (m.IsRoot()) {
      // root has concatenation of all colors and values
      if (colors.size() != values.size()) {
        throw std::runtime_error(
            "color-reduce: colors.size() != values.size()");
      }

      std::map<Scal, std::vector<Scal>> clmp;
      // reduce to map
      for (size_t k = 0; k < colors.size(); ++k) {
        auto cl = colors[k];
        auto& v = values[k];
        auto& vm = clmp[cl];
        vm.resize(v.size(), 0.);
        for (size_t i = 0; i < v.size(); ++i) {
          vm[i] += v[i];
        }
      }

      // divide by vf
      for (auto& it : clmp) {
        auto& v = it.second;
        Scal vf = v[0]; // XXX: assume vf is first
        Scal pi = M_PI;
        // XXX: assume r is second
        v[1] = std::pow(3. / (4. * pi) * vf, 1. / 3.);
        // divide remaining by vf
        for (size_t i = 2; i < v.size(); ++i) {
          v[i] /= vf;
        }
      }

      colors.clear();
      values.clear();
      for (auto& it : clmp) {
        colors.push_back(it.first);
        values.push_back(it.second);
      }
    }

    using TS = typename M::template OpCatT<Scal>;
    using TVS = typename M::template OpCatVT<Scal>;
    m.Bcast(std::make_shared<TS>(&colors));
    m.Bcast(std::make_shared<TVS>(&values));
  }
  if (sem("sphavg-init")) {
    if (var.Int["enable_shell"]) {
      sphavg.reset(new SA(m, var.Int["dim"]));
    }
  }
  if (sphavg && sem("sphavg-sph")) {
    const Scal shrr = var.Double["shell_rr"]; // shell inner radius relative
                                              // to equivalent radius
    const Scal shr = var.Double["shell_r"]; // shell inner radius absolute
    // shell total radius: rr * req + r
    const Scal shh = var.Double["shell_h"]; // shell thickness relative to h
    auto h = m.GetCellSize();

    vsph.clear();
    for (auto& s : values) {
      // XXX: adhoc, assume vf,r,x,y,z in values
      Vect x(s[2], s[3], s[4]);
      Scal r = s[1] * shrr + shr;
      vsph.emplace_back(x, r, h[0] * shh);
    }
  }
  if (sphavg && sem.Nested("sphavg-update")) {
    sphavg->Update(*fcvf[0], fcvel, fcvelm, dt, fcp, vsph);
  }
  if (sem("color-dump") && dm) {
    if (m.IsRoot()) {
      std::string s = GetDumpName("traj", ".csv", frame);
      std::cout << std::fixed << std::setprecision(8) << "dump"
                << " t=" << t << " to " << s << std::endl;
      std::ofstream o;
      o.open(s);
      o.precision(20);
      // header
      {
        o << "cl";
        for (size_t i = 0; i < names.size(); ++i) {
          o << "," << names[i];
        }
        o << std::endl;
      }
      // content
      for (size_t i = 0; i < colors.size(); ++i) {
        o << colors[i];
        for (auto v : values[i]) {
          o << "," << v;
        }
        o << "\n";
      }
    }
    if (sphavg) {
      if (m.IsRoot()) {
        std::string s = GetDumpName("trajsh", ".csv", frame);
        std::cout << std::fixed << std::setprecision(8) << "dump"
                  << " t=" << t << " to " << s << std::endl;
        std::ofstream o;
        o.open(s);
        o.precision(20);
        // header
        {
          auto nn = sphavg->GetNames();
          o << "cl";
          for (auto n : nn) {
            o << "," << n;
          }
          o << std::endl;
        }
        // content
        auto cl = colors;
        auto av = sphavg->GetAvg();
        if (cl.size() != av.size()) {
          throw std::runtime_error(
              "trajsh: cl.size()=" + std::to_string(cl.size()) +
              " != av.size()=" + std::to_string(av.size()));
        }
        for (size_t i = 0; i < cl.size(); ++i) {
          o << cl[i];
          for (auto& a : av[i].SerOut()) {
            o << "," << a;
          }
          o << "\n";
        }
      }
    }
  }
}

// ffst: force projections to append
// fcu: volume fraction
// fck: curvature
// ffsig: surface tension coefficient
template <class M>
void AppendSurfaceTension(
    const M& m, FieldFace<typename M::Scal>& ffst,
    const FieldCell<typename M::Scal>& fcu,
    const FieldCell<typename M::Scal>& fck,
    const FieldFace<typename M::Scal>& ffsig) {
  using Scal = typename M::Scal;
  for (auto f : m.Faces()) {
    const IdxCell cm = m.GetCell(f, 0);
    const IdxCell cp = m.GetCell(f, 1);
    const Scal um = fcu[cm];
    const Scal up = fcu[cp];
    const Scal k =
        (std::abs(um - 0.5) < std::abs(up - 0.5) ? fck[cm] : fck[cp]);
    const Scal hi = m.GetArea(f) / m.GetVolume(cp);
    const Scal ga = (up - um) * hi;
    if (ga != 0.) {
      ffst[f] += ga * (IsNan(k) ? 0 : k) * ffsig[f];
    }
  }
}

// ffst: force projections to append
// fcu: volume fraction
// fck: curvature
// ffsig: surface tension coefficient
template <class M>
void AppendSurfaceTension(
    const Embed<M>& eb, FieldFace<typename M::Scal>& ffst,
    const FieldCell<typename M::Scal>& fcu,
    const FieldCell<typename M::Scal>& fck,
    const FieldFace<typename M::Scal>& ffsig) {
  const auto& m = eb.GetMesh();
  using Scal = typename M::Scal;
  const auto fegp = eb.Gradient(fcu, MapCondFace(), 0, 0.);
  for (auto f : eb.Faces()) {
    const Scal gu = fegp[f];
    if (gu != 0.) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      const Scal um = fcu[cm];
      const Scal up = fcu[cp];
      const Scal k =
          (std::abs(um - 0.5) < std::abs(up - 0.5) ? fck[cm] : fck[cp]);
      ffst[f] += gu * (IsNan(k) ? 0 : k) * ffsig[f];
    }
  }
}

template <class M>
void AppendSurfaceTension(
    const M& m, FieldFace<typename M::Scal>& ffst, const GRange<size_t>& layers,
    const Multi<const FieldCell<typename M::Scal>*> fcu,
    const Multi<const FieldCell<typename M::Scal>*> fccl,
    const Multi<const FieldCell<typename M::Scal>*> fck,
    const FieldFace<typename M::Scal>& ffsig) {
  using Scal = typename M::Scal;
  constexpr Scal kClNone = -1;
  for (auto f : m.Faces()) {
    IdxCell cm = m.GetCell(f, 0);
    IdxCell cp = m.GetCell(f, 1);
    std::set<Scal> s;
    for (auto i : layers) {
      Scal clm = (*fccl[i])[cm];
      Scal clp = (*fccl[i])[cp];
      if (clm != kClNone) s.insert(clm);
      if (clp != kClNone) s.insert(clp);
    }
    for (auto cl : s) {
      Scal um = 0;
      Scal up = 0;
      Scal km = 0;
      Scal kp = 0;
      for (auto i : layers) {
        if ((*fccl[i])[cm] == cl) {
          um = (*fcu[i])[cm];
          km = (*fck[i])[cm];
        }
        if ((*fccl[i])[cp] == cl) {
          up = (*fcu[i])[cp];
          kp = (*fck[i])[cp];
        }
      }
      Scal k = (std::abs(um - 0.5) < std::abs(up - 0.5) ? km : kp);
      Scal hi = m.GetArea(f) / m.GetVolume(cp);
      Scal ga = (up - um) * hi;
      if (ga != 0.) {
        k = (IsNan(k) ? 0 : k);
        ffst[f] += ga * k * ffsig[f];
      }
    }
  }
}

template <class M>
void CalcSurfaceTension(
    const M& m, const GRange<size_t>& layers, const Vars& var,
    FieldCell<typename M::Vect>& fc_force,
    FieldFace<typename M::Scal>& ff_force,
    const FieldCell<typename M::Scal>& fc_sig, const MapCondFace& mf_sig,
    const Multi<const FieldCell<typename M::Scal>*> fck,
    const FieldCell<typename M::Scal>& fcvf,
    const FieldFace<typename M::Scal>& ffvfsm, const AdvectionSolver<M>* asb) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  auto st = var.String["surftens"];
  if (st == "div") { // divergence of tensor (Hu,Adam 2001)
    // volume fration gradient on cells
    const FieldCell<Vect> gc = Gradient(ffvfsm, m); // [s]
    // volume fration gradient on faces
    const FieldFace<Vect> gf =
        Interpolate(gc, GetCondZeroGrad<Vect>(mf_sig), m); // [i]
    auto stdiag = var.Double["stdiag"];
    for (auto c : m.Cells()) {
      Vect r(0);
      for (auto q : m.Nci(c)) {
        IdxFace f = m.GetFace(c, q);
        const auto& g = gf[f];
        // TODO: revise 1e-6
        auto n = g / (g.norm() + 1e-6); // inner normal
        auto s = m.GetOutwardSurface(c, q);
        r += s * (g.norm() * stdiag) - g * s.dot(n);
      }
      r /= m.GetVolume(c);
      // here: r = stdiag*div(|g|I) - div(g g/|g|)
      fc_force[c] += r * fc_sig[c];
    }
  } else if (st == "kn") { // curvature * normal
    FieldFace<Scal> ff_st(m, 0.); // surface tension projections
    const FieldFace<Scal> ff_sig = Interpolate(fc_sig, mf_sig, m);

    if (auto as = dynamic_cast<const Vofm<M>*>(asb)) {
      AppendSurfaceTension(
          m, ff_st, layers, as->GetFieldM(), as->GetColor(), fck, ff_sig);
    } else if (auto as = dynamic_cast<const VofEmbed<M>*>(asb)) {
      AppendSurfaceTension(
          as->GetEmbed(), ff_st, as->GetField(), *fck[0], ff_sig);
    } else if (auto as = dynamic_cast<const Vof<M>*>(asb)) {
      AppendSurfaceTension(m, ff_st, as->GetField(), *fck[0], ff_sig);
    } else {
      throw std::runtime_error("CalcSurfaceTension: unknown advection solver");
    }

    // zero on boundaries
    for (auto& it : mf_sig) {
      ff_st[it.first] = 0.;
    }

    // Surface tension decay between x0 and x1
    // XXX: adhoc TODO: revise
    const Scal x0 = var.Double["zerostx0"];
    const Scal x1 = var.Double["zerostx1"];
    // apply
    for (auto f : m.Faces()) {
      Scal x = m.GetCenter(f)[0];
      if (x > x0) {
        ff_st[f] *= std::max(0., (x1 - x) / (x1 - x0));
      }
    }

    // Append to force
    for (auto f : m.Faces()) {
      ff_force[f] += ff_st[f];
    }

    // Append Marangoni stress
    if (var.Int["marangoni"]) {
      if (auto as = dynamic_cast<const Vofm<M>*>(asb)) {
        using R = Reconst<Scal>;
        const auto fc_gsig = Gradient(ff_sig, m);
        const auto& fcn = *as->GetNormal()[0];
        const auto& fca = *as->GetAlpha()[0];
        const Vect h = m.GetCellSize();
        for (auto c : m.Cells()) {
          if (fcvf[c] > 0. && fcvf[c] < 1. && !IsNan(fca[c])) {
            Vect g = fc_gsig[c]; // sigma gradient
            Vect n = fcn[c] / fcn[c].norm(); // unit normal to interface
            Vect gt = g - n * g.dot(n);
            auto xx = R::GetCutPoly2(fcn[c], fca[c], h);
            Scal ar = std::abs(R::GetArea(xx, fcn[c]));
            Scal vol = h.prod();
            fc_force[c] += gt * (ar / vol);
          }
        }
      }
    }
  } else {
    throw std::runtime_error("Unknown surftens=" + st);
  }
}

template <class M>
void ProjectVolumeFlux(
    FieldFace<typename M::Scal>& ffv, const MapCondFaceFluid& mfc, M& m) {
  using Scal = typename M::Scal;
  using ExprFace = GVect<Scal, 3>;
  using Expr = GVect<Scal, M::dim * 2 + 2>;

  auto sem = m.GetSem();
  struct {
    FieldFace<ExprFace> ffe; // expression for corrected volume flux [i]
    FieldCell<Expr> fce; // linear system for pressure [i]
    FieldFace<bool> ffbd; // true for faces with boundary conditions
    FieldCell<Scal> fcp; // pressure (up to a constant)
  } * ctx(sem);
  auto& ffe = ctx->ffe;
  auto& fce = ctx->fce;
  auto& ffbd = ctx->ffbd;
  auto& fcp = ctx->fcp;

  if (sem("init")) {
    ffbd.Reinit(m, false);
    for (auto& p : mfc) {
      ffbd[p.first] = true;
    }

    ffe.Reinit(m);
    for (auto f : m.Faces()) {
      auto& e = ffe[f];
      if (!ffbd[f]) { // inner
        Scal a = -1; // XXX adhoc uniform
        e[0] = -a;
        e[1] = a;
      } else { // boundary
        e[0] = 0;
        e[1] = 0;
      }
      e[2] = ffv[f];
    }

    fce.Reinit(m, Expr(0));
    for (auto c : m.Cells()) {
      auto& e = fce[c];
      for (auto q : m.Nci(c)) {
        const IdxFace f = m.GetFace(c, q);
        const ExprFace v = ffe[f] * m.GetOutwardFactor(c, q);
        e[0] += v[1 - q % 2];
        e[1 + q] += v[q % 2];
        e[Expr::dim - 1] += v[2];
      }
    }
  }
  if (sem("solve")) {
    std::vector<Scal>* lsa;
    std::vector<Scal>* lsb;
    std::vector<Scal>* lsx;
    m.GetSolveTmp(lsa, lsb, lsx);
    lsx->resize(m.GetInBlockCells().size());
    size_t i = 0;
    for (auto c : m.Cells()) {
      (void)c;
      (*lsx)[i++] = 0;
    }
    auto l = ConvertLsCompact(fce, *lsa, *lsb, *lsx, m);
    using T = typename M::LS::T;
    l.t = T::symm;
    m.Solve(l);
  }
  if (sem("copy")) {
    std::vector<Scal>* lsa;
    std::vector<Scal>* lsb;
    std::vector<Scal>* lsx;
    m.GetSolveTmp(lsa, lsb, lsx);

    fcp.Reinit(m);
    size_t i = 0;
    for (auto c : m.Cells()) {
      fcp[c] = (*lsx)[i++];
    }
    m.Comm(&fcp);
  }
  if (sem("apply")) {
    for (auto f : m.Faces()) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      const auto& e = ffe[f];
      ffv[f] = e[0] * fcp[cm] + e[1] * fcp[cp] + e[2];
    }
  }
}

template <class M>
std::map<typename M::Scal, typename M::Scal> CalcArea(
    const GRange<size_t>& layers,
    const Multi<const FieldCell<typename M::Vect>*> fcn,
    const Multi<const FieldCell<typename M::Scal>*> fca,
    const Multi<const FieldCell<typename M::Scal>*> fccl, M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr Scal kClNone = -1; // TODO define kClNone once
  auto sem = m.GetSem();
  struct {
    std::vector<Scal> vcl; // color
    std::vector<Scal> vs; // area
  } * ctx(sem);
  auto& vcl = ctx->vcl;
  auto& vs = ctx->vs;
  if (sem("area")) {
    std::map<Scal, Scal> map;
    for (auto c : m.Cells()) {
      for (auto l : layers) {
        const Scal cl = (*fccl[l])[c];
        const Scal a = (*fca[l])[c];
        const Vect n = (*fcn[l])[c];
        if (cl != kClNone && !IsNan(a)) {
          using R = Reconst<Scal>;
          auto xx = R::GetCutPoly2(n, a, m.GetCellSize());
          map[cl] += std::abs(R::GetArea(xx, n));
        }
      }
    }
    for (auto p : map) {
      vcl.push_back(p.first);
      vs.push_back(p.second);
    }
    using T = typename M::template OpCatT<Scal>;
    m.Reduce(std::make_shared<T>(&vcl));
    m.Reduce(std::make_shared<T>(&vs));
  }
  if (sem("bcast")) {
    if (m.IsRoot()) {
      std::map<Scal, Scal> map;
      for (size_t i = 0; i < vcl.size(); ++i) {
        map[vcl[i]] += vs[i];
      }
      vcl.clear();
      vs.clear();
      for (auto p : map) {
        vcl.push_back(p.first);
        vs.push_back(p.second);
      }
    }
    using T = typename M::template OpCatT<Scal>;
    m.Bcast(std::make_shared<T>(&vcl));
    m.Bcast(std::make_shared<T>(&vs));
  }
  if (sem("map")) {
    std::map<Scal, Scal> map;
    for (size_t i = 0; i < vcl.size(); ++i) {
      map[vcl[i]] = vs[i];
    }
    return map;
  }
  return {};
}

template <class M>
std::map<typename M::Scal, typename M::Scal> CalcVolume(
    const GRange<size_t>& layers,
    const Multi<const FieldCell<typename M::Scal>*> fcu,
    const Multi<const FieldCell<typename M::Scal>*> fccl, M& m) {
  using Scal = typename M::Scal;
  static constexpr Scal kClNone = -1; // TODO define kClNone once
  auto sem = m.GetSem();
  struct {
    std::vector<Scal> vcl; // color
    std::vector<Scal> vvol; // volume
  } * ctx(sem);
  auto& vcl = ctx->vcl;
  auto& vvol = ctx->vvol;
  if (sem("volume")) {
    std::map<Scal, Scal> map;
    for (auto c : m.Cells()) {
      for (auto l : layers) {
        const Scal cl = (*fccl[l])[c];
        if (cl != kClNone) {
          map[cl] += (*fcu[l])[c] * m.GetVolume(c);
        }
      }
    }
    for (auto p : map) {
      vcl.push_back(p.first);
      vvol.push_back(p.second);
    }
    using T = typename M::template OpCatT<Scal>;
    m.Reduce(std::make_shared<T>(&vcl));
    m.Reduce(std::make_shared<T>(&vvol));
  }
  if (sem("bcast")) {
    if (m.IsRoot()) {
      std::map<Scal, Scal> map;
      for (size_t i = 0; i < vcl.size(); ++i) {
        map[vcl[i]] += vvol[i];
      }
      vcl.clear();
      vvol.clear();
      for (auto p : map) {
        vcl.push_back(p.first);
        vvol.push_back(p.second);
      }
    }
    using T = typename M::template OpCatT<Scal>;
    m.Bcast(std::make_shared<T>(&vcl));
    m.Bcast(std::make_shared<T>(&vvol));
  }
  if (sem("map")) {
    std::map<Scal, Scal> map;
    for (size_t i = 0; i < vcl.size(); ++i) {
      map[vcl[i]] = vvol[i];
    }
    return map;
  }
  return {};
}
