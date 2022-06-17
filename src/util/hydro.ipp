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
#include "func/init_bc.h"
#include "func/init_u.h"
#include "func/init_vel.h"
#include "func/primlist.h"
#include "parse/util.h"
#include "solver/approx.h"
#include "solver/approx_eb.h"
#include "solver/cond.h"
#include "solver/pois.h"
#include "solver/sphavg.h"
#include "solver/vof.h"
#include "solver/vofm.h"

#include "hydro.h"

template <class M>
void InitVel(FieldCell<typename M::Vect>& fcv, const Vars& var, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Vect2 = generic::Vect<Scal, 2>;
  const std::string vi = var.String["vel_init"];
  fcv.Reinit(m);
  if (auto ptr = ModuleInitVelocity<M>::GetInstance(vi)) {
    (*ptr)(fcv, var, m);
  } else if (vi == "grid_sin") {
    auto dx = var.Double["initvel_grid_sin_wavelength"];
    const Vect h = m.GetCellSize();
    for (auto c : m.AllCells()) {
      auto& v = fcv[c];
      auto x = m.GetCenter(c);
      v = Vect(0);
      v[1] = std::sin(x[0] / (dx * h[0]));
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
      auto z = (M::dim > 2 ? xx[2] : 0);
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
      auto z = (M::dim > 2 ? xx[2] : 0);
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
    Vect qc(0); // center
    qc[0] = (r1 + r0) * 0.5;
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
      Vect q(0);
      q[0] = xt;
      q[1] = xn;
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
      v[n.abs().argmin()] = 1;
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
      if (M::dim > 2) v[2] = vz.dot(x);
      fcv[c] = v;
    }
  } else if (vi == "pois" || vi == "poisy") {
    // Poiseuille with walls in y
    const Vect domain = m.GetGlobalLength();
    const Scal pv = var.Double["poisvel"]; // centerline velocity

    for (auto i : m.AllCells()) {
      const Scal y = m.GetCenter(i)[1] / domain[1];
      fcv[i][0] = y * (1. - y) * 4. * pv;
    }
  } else if (vi == "poisyz") {
    // Poiseuille with walls in y and z
    // Spiga 1994: Symmetric solution for velocity in rectangular ducts
    const Vect domain = m.GetGlobalLength();
    Scal mu = var.Double["poismu"]; // viscosity
    Scal pg = var.Double["poisgrad"]; // pressure gradient
    int im = var.Int["poisiter"]; // depth to evaluate series
    bool wym = var.Int["poiswym"]; // wallym
    bool wyp = var.Int["poiswyp"]; // wallyp
    bool wzm = var.Int["poiswzm"]; // wallzm
    bool wzp = var.Int["poiswzp"]; // wallzp
    Scal pi = M_PI;

    Scal ly = domain[1];
    Scal lz = M::dim > 2 ? domain[2] : 1;

    if ((!wym && !wyp) || (!wzm && !wzp)) {
      fassert(false, "poisyz: can't remove both walls");
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

    // TODO: tests for domain[1] != domain[2]
    for (auto i : m.AllCells()) {
      Scal y = m.GetCenter(i)[1] / ly;
      y = y + oy;
      Scal z = (M::dim > 2 ? m.GetCenter(i)[2] : 0) / lz;
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
      fcv[c] = Vect(Vect2(u, v)) * (y > h ? 0. : 1.);
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
      fcv[c] = Vect(Vect2(Ux, Uy));
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

      fcv[c] = Vect(Vect2(u, v));
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

      fcv[c] = Vect(Vect2(vx, vy)) * kvel;
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

      if (M::dim > 2) {
        fcv[c][2] = omz * omk;
      }
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

      fcv[c] = Vect(Vect2(u, v));
      if (y > h) {
        fcv[c] *= 0;
      }
    }
  } else if (vi == "list") {
    std::stringstream sstr(var.String["vellist_path"]);
    std::string fname;
    sstr >> fname;
    std::vector<typename UPrimList<Vect>::VelocityPrimitive> pp;

    if (fname == "inline") {
      if (m.IsRoot()) {
        std::cerr << __func__
                  << ": Reading inline list of primitives from vellist_path\n";
      }
      pp = UPrimList<Vect>::GetVelocityPrimitives(sstr, var.Int["dim"]);
    } else {
      if (m.IsRoot()) {
        std::cerr << __func__ << ": Open list of velocity primitives '" << fname
                  << "'\n";
      }
      std::ifstream fin(fname);
      fassert(fin.good(), "Can't open list of primitives '" + fname + "'");
      pp = UPrimList<Vect>::GetVelocityPrimitives(fin, var.Int["dim"]);
    }

    if (m.IsRoot()) {
      std::cerr << "Read " << pp.size() << " primitives" << std::endl;
    }

    fcv.Reinit(m, Vect(0));
    for (auto c : m.AllCells()) {
      Vect x = m.GetCenter(c);
      for (auto& p : pp) {
        fcv[c] += p.velocity(x);
      }
    }
  } else if (vi == "zero") {
    // nop
  } else {
    fassert(false, "Init(): unknown vel_init=" + vi);
  }
}

template <class Scal>
std::ostream& operator<<(std::ostream& out, const BCondAdvection<Scal>& fca) {
  using Halo = typename BCondAdvection<Scal>::Halo;
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
bool ParseAdvectionFaceCond(std::string str, BCondAdvection<Scal>& cfa) {
  using Halo = typename BCondAdvection<Scal>::Halo;
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

template <class MEB>
std::tuple<
    MapEmbed<BCondFluid<typename MEB::Vect>>,
    MapEmbed<BCondAdvection<typename MEB::Scal>>, MapEmbed<size_t>,
    std::vector<std::string>,
    std::vector<std::map<std::string, typename MEB::Scal>>>
InitBc(
    const Vars& var, const MEB& eb, std::set<std::string> known_keys,
    const FieldCell<bool>& fc_innermask) {
  using M = typename MEB::M;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using UI = UInitEmbedBc<M>;

  // default
  const Scal clear0 = var.Double["bcc_clear0"];
  const Scal clear1 = var.Double["bcc_clear1"];
  const Scal inletcl = var.Double["inletcl"];
  const Scal fill_vf = var.Double["bcc_fill"];

  auto default_adv = [clear0, clear1, inletcl,
                      fill_vf](const BCondFluid<Vect>& bcf) {
    using Halo = typename BCondAdvection<Scal>::Halo;
    static constexpr Scal kClNone = -1; // TODO define kClNone once
    BCondAdvection<Scal> ca;
    ca.nci = bcf.nci;
    switch (bcf.type) {
      case BCondFluidType::wall:
      case BCondFluidType::slipwall:
        ca.halo = Halo::fill;
        ca.clear0 = clear0;
        ca.clear1 = clear1;
        ca.fill_vf = fill_vf;
        ca.fill_cl = kClNone;
        break;
      case BCondFluidType::inlet:
      case BCondFluidType::inletflux:
      case BCondFluidType::inletpressure:
        ca.halo = Halo::fill;
        ca.clear0 = 0;
        ca.clear1 = 1;
        ca.fill_vf = 0;
        ca.fill_cl = inletcl;
        break;
      case BCondFluidType::outlet:
      case BCondFluidType::outletpressure:
        ca.halo = Halo::fill;
        ca.clear0 = 1;
        ca.clear1 = 1;
        ca.fill_vf = 0;
        ca.fill_cl = kClNone;
        break;
      case BCondFluidType::symm:
        ca.halo = Halo::reflect;
        break;
    }
    return ca;
  };

  auto parse_key_value = [&](std::string s) {
    std::stringstream buf(s);
    std::string key;
    Scal value;
    buf >> key;
    buf >> value;
    return std::pair<std::string, Scal>({key, value});
  };

  auto parse = [&](std::string list, size_t nci, Vect face_center,
                   Vect face_normal) {
    std::vector<std::string> ss_adv; // strings to try as advection conditions
    bool found_fluid = false;
    std::map<std::string, Scal> custom;
    BCondFluid<Vect> bc_fluid;
    for (std::string s : SplitByDelimiter(list, ',')) {
      auto key_value = parse_key_value(s);
      if (known_keys.count(key_value.first)) {
        custom.insert(key_value);
        continue;
      }
      auto p = ParseBCondFluid<Vect>(s, nci, face_center, face_normal);
      if (p.first) { // parsed as fluid condition
        found_fluid = true;
        bc_fluid = p.second;
      } else { // try as advection later
        ss_adv.push_back(s);
      }
    }
    fassert(found_fluid, "No fluid condition found in '" + list + "'");
    auto bc_adv = default_adv(bc_fluid);
    for (std::string s : ss_adv) {
      fassert(
          ParseAdvectionFaceCond(s, bc_adv),
          "No advection condition found in '" + s + "'");
    }
    return std::make_tuple(bc_fluid, bc_adv, custom);
  };

  MapEmbed<BCondFluid<Vect>> me_fluid;
  MapEmbed<BCondAdvection<Scal>> me_adv;

  std::stringstream buf;
  {
    std::stringstream path(var.String["bc_path"]);
    std::string filename;
    path >> filename;
    if (filename == "inline") {
      buf << path.rdbuf();
    } else {
      filename = path.str();
      std::ifstream fin(filename);
      fassert(fin.good(), "Can't open boundary conditions '" + filename + "'");
      buf << fin.rdbuf();
    }
  }
  std::vector<std::string> vdesc;
  MapEmbed<size_t> me_group;
  MapEmbed<size_t> me_nci;
  std::tie(me_group, me_nci, vdesc) = UI::ParseGroups(buf, eb, fc_innermask);
  std::vector<std::map<std::string, Scal>> vcustom(vdesc.size());

  me_group.LoopPairs([&](const auto& p) {
    const auto cf = p.first;
    const auto group = p.second;
    std::tie(me_fluid[cf], me_adv[cf], vcustom[group]) = parse(
        vdesc[group], me_nci.at(cf), eb.GetFaceCenter(cf), eb.GetNormal(cf));
  });
  return {me_fluid, me_adv, me_group, vdesc, vcustom};
}

template <class MEB>
void DumpBcPoly(
    const std::string filename, const MapEmbed<size_t>& me_group,
    const MapEmbed<typename MEB::Scal>& me_contang, const MEB& meb,
    typename MEB::M& m) {
  using M = typename MEB::M;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  auto sem = m.GetSem("dumppoly");
  struct {
    std::vector<std::vector<Vect>> dpoly;
    std::vector<Scal> dgroup;
    std::vector<Scal> dcontang;
    std::vector<Scal> dface;
  } * ctx(sem);
  auto& dpoly = ctx->dpoly;
  auto& dgroup = ctx->dgroup;
  auto& dcontang = ctx->dcontang;
  auto& dface = ctx->dface;
  if (sem("local")) {
    me_group.LoopPairs([&](auto p) {
      const auto cf = p.first;
      dpoly.push_back(meb.GetFacePoly(cf));
      dgroup.push_back(p.second);
      dcontang.push_back(me_contang.at(cf));
      dface.push_back(meb.IsCell(cf) ? 0 : 1);
    });
    m.Reduce(&dpoly, Reduction::concat);
    m.Reduce(&dgroup, Reduction::concat);
    m.Reduce(&dcontang, Reduction::concat);
    m.Reduce(&dface, Reduction::concat);
  }
  if (sem("write")) {
    if (m.IsRoot()) {
      dump::Vtk<Vect>::WriteVtkPoly(
          filename, dpoly, nullptr, {&dgroup, &dcontang, &dface},
          {"group", "contang", "face"}, "Boundary conditions", true, true,
          true);
    }
  }
}

template <class M>
void GetFluidCellCond(
    const Vars& var, M& m, MapCell<std::shared_ptr<CondCellFluid>>& mcvel) {
  using Vect = typename M::Vect;
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
      m.Reduce(&pdist, Reduction::minloc);
    }
  }

  if (sem("set")) {
    // Fix pressure at one cell
    if (auto* p = var.Double.Find("pfixed")) {
      if (pdist.second == m.GetId()) {
        Vect x(var.Vect["pfixed_x"]);
        IdxCell c = m.FindNearestCell(x);
        mcvel[c] = std::make_shared<fluid_condition::GivenPressureFixed<M>>(*p);
        std::cerr << "pfixed id=" << pdist.second << " dist=" << pdist.first
                  << std::endl;
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

// Computes velocity fcvel from vorticity fcvort
// Scalar vorticity is stored in fcvort[c][0].
// Component fcvort[c][1] is ignored.
template <class M>
void InitVort(
    const FieldCell<typename M::Vect>& fcvort,
    FieldCell<typename M::Vect>& fcvel, FieldCell<typename M::Vect>* fcpot,
    const MapEmbed<BCondFluid<typename M::Vect>>& mebc_fluid,
    std::shared_ptr<linear::Solver<M>> linsolver, M& m,
    generic::Vect<typename M::Scal, 2>*, bool zero_dirichlet) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  auto sem = m.GetSem("initvort");
  using UEB = UEmbed<M>;
  struct {
    MapEmbed<BCond<Vect>> mebc_vel; // velocity conditions
    MapEmbed<BCond<Scal>> mebc_pot; // potential conditions
    FieldCell<Scal> fcvort_scal; // one component of vorticity
    FieldCell<Scal> fcpot; // velocity potential
  } * ctx(sem);
  auto& t = *ctx;
  if (sem("initpois")) {
    t.fcpot.Reinit(m);
    t.mebc_vel = ConvertBCondFluidToVelocity<M>(mebc_fluid);
    t.mebc_pot = GetBCondZeroGrad<Scal>(t.mebc_vel);
    if (zero_dirichlet) {
      mebc_fluid.LoopPairs([&](auto p) {
        const auto cf = p.first;
        const auto bc = p.second;
        t.mebc_pot[cf] = BCond<Scal>(BCondType::dirichlet, bc.nci, 0);
      });
    }
  }
  if (sem("init")) {
    t.fcvort_scal = GetComponent(fcvort, 0);
    for (auto c : m.Cells()) {
      t.fcvort_scal[c] *= -1;
    }
  }
  if (sem.Nested("solve")) {
    SolvePoisson(t.fcpot, t.fcvort_scal, t.mebc_pot, linsolver, m, false);
  }
  if (sem("vel")) {
    auto ffg = UEB::Gradient(t.fcpot, t.mebc_pot, m);
    auto fcg = UEB::AverageGradient(ffg, m);
    for (auto c : m.Cells()) {
      fcvel[c][0] = fcg[c][1];
      fcvel[c][1] = -fcg[c][0];
    }
    if (fcpot) {
      fcpot->Reinit(m);
      SetComponent(*fcpot, 0, t.fcpot);
    }
    m.Comm(&fcvel);
  }
}

// Computes velocity fcvel from vorticity fcvort
template <class M>
void InitVort(
    const FieldCell<typename M::Vect>& fcvort,
    FieldCell<typename M::Vect>& fcvel, FieldCell<typename M::Vect>* fcpot,
    const MapEmbed<BCondFluid<typename M::Vect>>& mebc_fluid,
    std::shared_ptr<linear::Solver<M>> linsolver, M& m,
    generic::Vect<typename M::Scal, 3>*, bool zero_dirichlet) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  auto sem = m.GetSem("initvort");
  struct {
    MapEmbed<BCond<Vect>> mebc_vel; // velocity conditions
    MapEmbed<BCond<Scal>> mebc_pot_scal; // conditions for potential
    FieldCell<Scal> fcvort_scal; // one component of vorticity
    FieldCell<Vect> fcpot; // velocity potential
    FieldCell<Scal> fcpot_scal; // one component of potential
  } * ctx(sem);
  auto& t = *ctx;
  if (sem("initpois")) {
    t.fcpot.Reinit(m);
    t.mebc_vel = ConvertBCondFluidToVelocity<M>(mebc_fluid);
    t.mebc_pot_scal = GetBCondZeroGrad<Scal>(t.mebc_vel);
    if (zero_dirichlet) {
      mebc_fluid.LoopPairs([&](auto p) {
        const auto cf = p.first;
        const auto bc = p.second;
        t.mebc_pot_scal[cf] = BCond<Scal>(BCondType::dirichlet, bc.nci, 0);
      });
    }
  }
  for (size_t d = 0; d < M::dim; ++d) {
    const std::string dname = std::to_string(d);
    if (sem("init-" + dname)) {
      t.fcvort_scal = GetComponent(fcvort, d);
      for (auto c : m.Cells()) {
        t.fcvort_scal[c] *= -1;
      }
    }
    if (sem.Nested("solve-" + dname)) {
      SolvePoisson(t.fcpot_scal, t.fcvort_scal, t.mebc_pot_scal, linsolver, m);
    }
    if (sem("post-" + dname)) {
      SetComponent(t.fcpot, d, t.fcpot_scal);
    }
  }
  if (sem("vel")) {
    fcvel = UEmbed<M>::GetVort(t.fcpot, t.mebc_vel, m);
    if (fcpot) {
      fcpot->Reinit(m);
      SetComponent(*fcpot, 0, t.fcpot_scal);
    }
    m.Comm(&fcvel);
  }
}

template <class M>
void InitVort(
    const FieldCell<typename M::Vect>&, FieldCell<typename M::Vect>&,
    FieldCell<typename M::Vect>*, const MapEmbed<BCondFluid<typename M::Vect>>&,
    std::shared_ptr<linear::Solver<M>>, M&, generic::Vect<typename M::Scal, 4>*,
    bool) {}

// Computes velocity fcvel from vorticity fcvort
template <class M>
void InitVort(
    const FieldCell<typename M::Vect>& fcvort,
    FieldCell<typename M::Vect>& fcvel, FieldCell<typename M::Vect>* fcpot,
    const MapEmbed<BCondFluid<typename M::Vect>>& mebc_fluid,
    std::shared_ptr<linear::Solver<M>> linsolver, M& m, bool zero_dirichlet) {
  InitVort(
      fcvort, fcvel, fcpot, mebc_fluid, linsolver, m,
      (typename M::Vect*)nullptr, zero_dirichlet);
}

template <class EB>
void CalcTraj(
    EB& eb, const GRange<size_t>& layers,
    const Multi<const FieldCell<typename EB::Scal>*>& fcvf,
    const Multi<const FieldCell<typename EB::Scal>*>& fccl,
    const Multi<const FieldCell<typename EB::MIdx>*>& fcim,
    const FieldCell<typename EB::Scal>& fcp,
    const FieldCell<typename EB::Vect>& fcvel,
    /*out*/
    std::vector<std::string>& column_names,
    std::vector<typename EB::Scal>& row_colors,
    std::vector<std::vector<typename EB::Scal>>& table) {
  using M = typename EB::M;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  constexpr Scal kClNone = -1;
  auto& m = eb.GetMesh();

  auto sem = m.GetSem("calctraj");
  if (sem("color-calc")) {
    std::map<Scal, std::vector<Scal>> cl2row;

    // add scalar name
    auto add_name = [&](const std::string nm) { column_names.push_back(nm); };
    // add vector name
    auto add_name_vect = [&](const std::string nm) {
      column_names.push_back(nm + "x");
      column_names.push_back(nm + "y");
      column_names.push_back(nm + "z");
    };

    // XXX: adhoc, the following order assumed in post: vf,r,x,y,z,...

    // list of vars
    column_names.clear();
    add_name("vf");
    add_name("r");
    add_name_vect("");
    add_name_vect("v");
    add_name("p");
    add_name("xx");
    add_name("xy");
    add_name("xz");
    add_name("yy");
    add_name("yz");
    add_name("zz");

    // traverse cells, append to cl2row
    for (auto c : eb.Cells()) {
      for (auto l : layers) {
        auto cl = (*fccl[l])[c];
        if (cl != kClNone) {
          auto& v = cl2row[cl]; // vector for data
          auto x = m.GetCenter(c); // cell center
          // translation by image vector
          x += Vect((*fcim[l])[c]) * m.GetGlobalLength();
          const auto w = (*fcvf[l])[c] * m.GetVolume(c); // volume

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
            for (auto d : M::dirs) {
              add(a[d]);
            }
          };

          // list of vars, XXX: keep consistent with column_names
          add(w); // vf, must be first, divided  after reduction
          add(0); // r, must be second, computed after reduction
          addv(x * w); // x
          addv(fcvel[c] * w); // velocity
          add(fcp[c] * w); // pressure
          add(x[0] * x[0] * w); // xx
          add(x[0] * x[1] * w); // xy
          if (M::dim > 2) {
            add(x[0] * x[2] * w); // xz
          }
          add(x[1] * x[1] * w); // yy
          if (M::dim > 2) {
            add(x[1] * x[2] * w); // yz
            add(x[2] * x[2] * w); // zz
          }
        }
      }
    }
    // copy to vector
    row_colors.clear();
    table.clear();
    for (auto& it : cl2row) {
      row_colors.push_back(it.first);
      table.push_back(it.second);
    }
    m.Reduce(&row_colors, Reduction::concat);
    m.Reduce(&table, Reduction::concat);
  }
  if (sem("color-post")) {
    if (m.IsRoot()) {
      // root has concatenation of all row_colors and table
      fassert_equal(row_colors.size(), table.size());

      std::map<Scal, std::vector<Scal>> cl2row;
      // reduce to map
      for (size_t k = 0; k < row_colors.size(); ++k) {
        auto cl = row_colors[k];
        auto& v = table[k];
        auto& vm = cl2row[cl];
        vm.resize(v.size(), 0.);
        for (size_t i = 0; i < v.size(); ++i) {
          vm[i] += v[i];
        }
      }

      // divide by vf
      for (auto& it : cl2row) {
        auto& v = it.second;
        Scal vf = v[0]; // XXX: assume vf is first
        const Scal pi = M_PI;
        // XXX: assume r is second
        if (m.GetEdim() == 3) {
          v[1] = std::pow(3. / (4. * pi) * vf, 1. / 3.);
        } else { // dim 2
          v[1] = std::sqrt(vf / pi / (M::dim > 2 ? m.GetCellSize()[2] : 1));
        }
        // divide remaining by vf
        for (size_t i = 2; i < v.size(); ++i) {
          v[i] /= vf;
        }
      }

      row_colors.clear();
      table.clear();
      for (auto& it : cl2row) {
        row_colors.push_back(it.first);
        table.push_back(it.second);
      }
    }

    m.Bcast(&row_colors);
    m.Bcast(&table);
  }
}

template <class EB>
void DumpTraj(
    EB& eb, bool dm, const Vars& var, size_t frame, typename EB::Scal time,
    const GRange<size_t>& layers,
    const Multi<const FieldCell<typename EB::Scal>*>& fcvf,
    const Multi<const FieldCell<typename EB::Scal>*>& fccl,
    const Multi<const FieldCell<typename EB::MIdx>*>& fcim,
    const FieldCell<typename EB::Scal>& fcp,
    const FieldCell<typename EB::Vect>& fcvel,
    const FieldCell<typename EB::Vect>& fcvelm, typename EB::Scal dt) {
  using M = typename EB::M;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Vect3 = generic::Vect<Scal, 3>;
  using SA = Sphavg<M>; // spherical averages
  auto& m = eb.GetMesh();

  auto sem = m.GetSem("dumptraj");
  struct {
    std::vector<std::string> column_names; // color reduce: variable name
    std::vector<Scal> row_colors; // color reduce: cl
    std::vector<std::vector<Scal>> table; // color reduce: vector
    std::unique_ptr<SA> sphavg; // spherical averages
    std::vector<typename SA::Sph> vsph;
  } * ctx(sem);
  auto& t = *ctx;
  auto& sphavg = ctx->sphavg;
  auto& vsph = ctx->vsph;
  if (sem.Nested("calc")) {
    CalcTraj(
        eb, layers, fcvf, fccl, fcim, fcp, fcvel, t.column_names, t.row_colors,
        t.table);
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
    for (auto& s : t.table) {
      // XXX: adhoc, assume vf,r,x,y,z in t.table
      const Vect x(Vect3(s[2], s[3], s[4]));
      const Scal r = s[1] * shrr + shr;
      vsph.emplace_back(x, r, h[0] * shh);
    }
  }
  if (sphavg && sem.Nested("sphavg-update")) {
    sphavg->Update(*fcvf[0], fcvel, fcvelm, dt, fcp, vsph);
  }
  if (sem("color-dump") && dm) {
    if (m.IsRoot()) {
      std::string s = GetDumpName("traj", ".csv", frame);
      std::cerr << std::fixed << std::setprecision(8) << "dump"
                << " t=" << time << " to " << s << std::endl;
      std::ofstream o;
      o.open(s);
      o.precision(20);
      // header
      {
        o << "cl";
        for (size_t i = 0; i < t.column_names.size(); ++i) {
          o << "," << t.column_names[i];
        }
        o << std::endl;
      }
      // content
      for (size_t i = 0; i < t.row_colors.size(); ++i) {
        o << t.row_colors[i];
        for (auto v : t.table[i]) {
          o << "," << v;
        }
        o << "\n";
      }
    }
    if (sphavg) {
      if (m.IsRoot()) {
        std::string s = GetDumpName("trajsh", ".csv", frame);
        std::cerr << std::fixed << std::setprecision(8) << "dump"
                  << " t=" << time << " to " << s << std::endl;
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
        const auto& avg = sphavg->GetAvg();
        fassert_equal(t.row_colors.size(), avg.size());
        for (size_t i = 0; i < t.row_colors.size(); ++i) {
          o << t.row_colors[i];
          for (auto& a : avg[i].SerOut()) {
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
    const M& eb, FieldFace<typename M::Scal>& ffst,
    const FieldCell<typename M::Scal>& fcu,
    const FieldCell<typename M::Scal>& fck,
    const FieldFace<typename M::Scal>& ffsig) {
  using Scal = typename M::Scal;
  const Scal h = eb.GetCellSize()[0];
  for (auto f : eb.Faces()) {
    const IdxCell cm = eb.GetCell(f, 0);
    const IdxCell cp = eb.GetCell(f, 1);
    const Scal um =
        std::max(0., std::min(1., fcu[cm] / eb.GetVolumeFraction(cm)));
    const Scal up =
        std::max(0., std::min(1., fcu[cp] / eb.GetVolumeFraction(cp)));
    const Scal ga = (up - um) / h;
    if (ga != 0.) {
      Scal k = (std::abs(um - 0.5) < std::abs(up - 0.5) ? fck[cm] : fck[cp]);
      k = (IsNan(k) ? 0 : k);
      ffst[f] += ga * k * ffsig[f];
    }
  }
}

template <class M>
void AppendSurfaceTension(
    const M& eb, FieldFace<typename M::Scal>& ffst,
    const GRange<size_t>& layers,
    const Multi<const FieldCell<typename M::Scal>*> fcu,
    const Multi<const FieldCell<typename M::Scal>*> fccl,
    const Multi<const FieldCell<typename M::Scal>*> fck,
    const FieldFace<typename M::Scal>& ffsig) {
  using Scal = typename M::Scal;
  constexpr Scal kClNone = -1;
  const Scal h = eb.GetCellSize()[0];
  for (auto f : eb.Faces()) {
    const IdxCell cm = eb.GetCell(f, 0);
    const IdxCell cp = eb.GetCell(f, 1);
    std::set<Scal> s;
    for (auto i : layers) {
      const Scal clm = (*fccl[i])[cm];
      const Scal clp = (*fccl[i])[cp];
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
      um = std::max(0., std::min(1., um / eb.GetVolumeFraction(cm)));
      up = std::max(0., std::min(1., up / eb.GetVolumeFraction(cp)));
      const Scal ga = (up - um) / h;
      if (ga != 0.) {
        Scal k = (std::abs(um - 0.5) < std::abs(up - 0.5) ? km : kp);
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
    const FieldCell<typename M::Scal>& fc_sig,
    const MapEmbed<BCond<typename M::Scal>>& mf_sig,
    const Multi<const FieldCell<typename M::Scal>*> fck,
    const FieldCell<typename M::Scal>& fcvf,
    const FieldFace<typename M::Scal>& ffvfsm, const AdvectionSolver<M>* asb) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  auto st = var.String["surftens"];
  if (st == "div") { // divergence of tensor (Hu,Adam 2001)
    // volume fration gradient on cells
    const FieldCell<Vect> gc = UEmbed<M>::AverageGradient(ffvfsm, m); // [s]
    // volume fration gradient on faces
    const FieldFace<Vect> gf =
        UEmbed<M>::Interpolate(gc, GetBCondZeroGrad<Vect>(mf_sig), m);
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
    const FieldFace<Scal> ff_sig = UEmbed<M>::Interpolate(fc_sig, mf_sig, m);

    bool found = false; // separate if's needed to prevent Wshadow
    if (auto* as = dynamic_cast<const Vofm<M>*>(asb)) {
      found = true;
      AppendSurfaceTension(
          m, ff_st, layers, as->GetFieldM(), as->GetColor(), fck, ff_sig);
    }
    if (auto* as = dynamic_cast<const Vofm<Embed<M>>*>(asb)) {
      found = true;
      AppendSurfaceTension(
          as->GetEmbed(), ff_st, layers, as->GetFieldM(), as->GetColor(), fck,
          ff_sig);
    }
    if (auto* as = dynamic_cast<const Vof<M>*>(asb)) {
      found = true;
      AppendSurfaceTension(m, ff_st, as->GetField(), *fck[0], ff_sig);
    }
    if (auto* as = dynamic_cast<const Vof<Embed<M>>*>(asb)) {
      found = true;
      AppendSurfaceTension(
          as->GetEmbed(), ff_st, as->GetField(), *fck[0], ff_sig);
    }
    if (!found) {
      fassert(false, "CalcSurfaceTension: unknown advection solver");
    }

    // zero on boundaries
    for (auto& it : mf_sig.GetMapFace()) {
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
        const auto fc_gsig = UEmbed<M>::AverageGradient(ff_sig, m);
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
    fassert(false, "Unknown surftens=" + st);
  }
}

template <class M>
void ProjectVolumeFlux(
    FieldFace<typename M::Scal>& ffv,
    const MapEmbed<BCondFluid<typename M::Vect>>& mfc,
    std::shared_ptr<linear::Solver<M>> linsolver, M& m) {
  using Scal = typename M::Scal;
  using ExprFace = generic::Vect<Scal, 3>;
  using Expr = generic::Vect<Scal, M::dim * 2 + 2>;

  auto sem = m.GetSem();
  struct {
    FieldFace<ExprFace> ffe; // expression for corrected volume flux [i]
    FieldCell<Expr> fc_system; // linear system for pressure [i]
    FieldFace<bool> ffbd; // true for faces with boundary conditions
    FieldCell<Scal> fcp; // pressure (up to a constant)
  } * ctx(sem);
  auto& ffe = ctx->ffe;
  auto& fc_system = ctx->fc_system;
  auto& ffbd = ctx->ffbd;
  auto& fcp = ctx->fcp;
  auto& t = *ctx;

  if (sem("init")) {
    ffbd.Reinit(m, false);
    for (auto& p : mfc.GetMapFace()) {
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

    fc_system.Reinit(m, Expr(0));
    for (auto c : m.Cells()) {
      auto& e = fc_system[c];
      for (auto q : m.Nci(c)) {
        const IdxFace f = m.GetFace(c, q);
        m.AppendExpr(e, ffe[f] * m.GetOutwardFactor(c, q), q);
      }
    }
  }
  if (sem.Nested("solve")) {
    linsolver->Solve(t.fc_system, nullptr, t.fcp, m);
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
    m.Reduce(&vcl, Reduction::concat);
    m.Reduce(&vs, Reduction::concat);
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
    m.Bcast(&vcl);
    m.Bcast(&vs);
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
    m.Reduce(&vcl, Reduction::concat);
    m.Reduce(&vvol, Reduction::concat);
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
    m.Bcast(&vcl);
    m.Bcast(&vvol);
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
