// Created by Petr Karnakov on 11.07.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>
#include <numeric>

#include "geom/vect.h"

template <class Scal_>
class Reconst {
 public:
  using Scal = Scal_;
  using Vect2 = generic::Vect<Scal, 2>;
  using Vect3 = generic::Vect<Scal, 3>;
  using Vect4 = generic::Vect<Scal, 4>;

  static Scal Clip(Scal a, Scal min, Scal max) {
    return std::max(min, std::min(max, a));
  }

  static Scal cube(Scal a) {
    return a * a * a;
  }

  static Scal sqr(Scal a) {
    return a * a;
  }

  // Returns k-th root of equation `ax^3 + bx^2 + cx + d = 0`
  static Scal SolveCubic(Scal a, Scal b, Scal c, Scal d, int k) {
    Scal p = (3. * a * c - b * b) / (3. * a * a);
    p = std::min(p, 0.);
    Scal q =
        (2. * cube(b) - 9. * a * b * c + 27. * a * a * d) / (27. * cube(a));
    Scal r = 3. * q * std::sqrt(-3. / p) / (2. * p);
    r = std::max(-1., std::min(1., r));
    Scal t = 2. * std::sqrt(-p / 3.) *
             std::cos(1. / 3. * std::acos(r) - 2. * M_PI * k / 3.);
    Scal x = t - b / (3. * a);
    return x;
  }

  // Returns center of a cloud of points
  template <class Vect>
  static Vect GetCenter(const std::vector<Vect>& xx) {
    Vect res(0);
    for (auto& x : xx) {
      res += x;
    }
    return res / xx.size();
  }

  // Nearest point to line between ends.
  // x: target point
  // x0,x1: line ends
  template <class Vect>
  static Vect GetNearest(Vect x, Vect x0, Vect x1) {
    const Vect l = x1 - x0;
    Scal k = l.dot(x - x0) / l.sqrnorm();
    return x0 + l * Clip(k, 0, 1);
  }

  static void Sort(Vect2& v) {
    if (v[1] < v[0]) {
      std::swap(v[0], v[1]);
    }
  }
  static void Sort(Vect3& v) {
    if (v[1] < v[0]) {
      std::swap(v[0], v[1]);
    }
    if (v[2] < v[1]) {
      std::swap(v[1], v[2]);
    }
    if (v[1] < v[0]) {
      std::swap(v[0], v[1]);
    }
  }
  static void Sort(Vect4& v) {
    std::sort(&v[0], &v[v.size()]);
  }

  // Returns sequence of indices `r` such that v[r] is sorted
  static generic::Vect<size_t, 2> Argsort(Vect2 v) {
    generic::Vect<size_t, 2> r(0ul, 1ul);
    if (v[1] < v[0]) {
      std::swap(v[0], v[1]);
      std::swap(r[0], r[1]);
    }
    return r;
  }
  static generic::Vect<size_t, 3> Argsort(Vect3 v) {
    generic::Vect<size_t, 3> r(0ul, 1ul, 2ul);
    if (v[1] < v[0]) {
      std::swap(v[0], v[1]);
      std::swap(r[0], r[1]);
    }
    if (v[2] < v[1]) {
      std::swap(v[1], v[2]);
      std::swap(r[1], r[2]);
    }
    if (v[1] < v[0]) {
      std::swap(v[0], v[1]);
      std::swap(r[0], r[1]);
    }
    return r;
  }
  static generic::Vect<size_t, 4> Argsort(Vect4 v) {
    std::array<size_t, 4> r;
    std::iota(r.begin(), r.end(), 0);
    std::stable_sort(r.begin(), r.end(), [begin = &v[0]](auto a, auto b) {
      return *(begin + a) < *(begin + b);
    });
    return generic::Vect<size_t, 4>(&r[0]);
  }

  // Volume fraction from plane constant in unit cell with ordered normal.
  // nx,ny: normal, 0 <= nx <= ny
  // a: plane constant, -0.5*(nx + ny) <= a <= 0
  static Scal GetLineU0(Vect2 n, Scal a) {
    const Scal nx = n[0];
    const Scal ny = n[1];
    const Scal f = a + 0.5 * (nx + ny);
    if (nx > f) {
      return sqr(f) / (2 * nx * ny);
    }
    return a / ny + 0.5;
  }
  // Volume fraction from plane constant in unit cell with ordered normal.
  // nx,ny,nz: normal, 0 <= nx <= ny <= nz
  // a: plane constant, -0.5*n.sum() <= a <= 0
  static Scal GetLineU0(Vect3 n, Scal a) {
    const Scal nx = n[0];
    const Scal ny = n[1];
    const Scal nz = n[2];
    Scal f = a + 0.5 * (nx + ny + nz);

    if (f <= 0) {
      return 0;
    }
    if (nx >= f) { // f>0, nx>0
      return cube(f) / (6 * nx * ny * nz);
    }
    if (ny >= f) {
      return (3 * sqr(f) - 3 * f * nx + sqr(nx)) / (6 * ny * nz);
    }
    if (nz >= f) {
      if (nx + ny >= f) { // ny>=f/2, ny<f, nx>=f-ny, nx>0
        return (3 * sqr(f) - 3 * f * nx + sqr(nx) -
                std::min(1., (f - ny) / nx) * sqr(f - ny)) /
               (6 * ny * nz);
      }
      return (2 * f - nx - ny) / (2 * nz);
    }
    return (cube(f) - cube(f - nx) - cube(f - ny) - cube(f - nz)) /
           (6 * nx * ny * nz);
  }
  // Volume fraction from plane constant in unit cell with ordered normal.
  // nx,ny,nz: normal, 0 <= nx <= ny <= nz
  // a: plane constant, -0.5*n.sum() <= a <= 0
  static Scal GetLineU0(Vect4 n, Scal a) {
    (void)n;
    (void)a;
    return 0; // XXX not implemented
  }

  // Volume fraction from plane constant in unit cell.
  // n0: normal
  // a: plane constant
  // Returns:
  // u: volume fraction
  // Equation of reconstructed plane
  // x.dot(n) = a
  template <class Vect>
  static Scal GetLineU1(Vect n0, Scal a) {
    Vect n = n0.abs();
    Sort(n);
    a = Clip(a, -0.5 * n.sum(), 0.5 * n.sum());

    if (a < 0) {
      return GetLineU0(n, a);
    }
    return 1 - GetLineU0(n, -a);
  }

  // Volume fraction from line constant in rectangular cell.
  // n : normal
  // a: line constant
  // h: cell size
  // Returns:
  // u: volume fraction
  template <class Vect>
  static Scal GetLineU(Vect n, Scal a, Vect h) {
    return GetLineU1(n * h, a);
  }

  // Plane constant from volume fraction in unit cell.
  // nx,ny,nz: normal, 0 <= nx <= ny <= nz
  // u: volume fraction, 0 <= u <= 0.5
  // Returns:
  // a: plane constant
  static Scal GetLineA0(Vect2 n, Scal u) {
    const Scal nx = n[0];
    const Scal ny = n[1];
    if (2 * ny * u <= nx) {
      return std::sqrt(2 * nx * ny * u) - 0.5 * (nx + ny);
    }
    return ny * (u - 0.5);
  }
  static Scal GetLineA0(Vect3 n, Scal u) {
    const Scal nx = n[0];
    const Scal ny = n[1];
    const Scal nz = n[2];
    Scal f;
    if (6. * ny * nz * u < sqr(nx)) {
      f = std::pow(6. * nx * ny * nz * u, 1. / 3.);
    } else if (6. * ny * nz * u < 3. * sqr(ny) - 3 * nx * ny + sqr(nx)) {
      f = 0.5 * nx + std::sqrt(2. * ny * nz * u - sqr(nx) / 12.);
    } else if (nz > nx + ny) {
      if (2. * nz * u < nx + ny) {
        f = SolveCubic(
            1., -3. * (nx + ny), 3. * (sqr(nx) + sqr(ny)),
            -(cube(nx) + cube(ny)) + 6. * nx * ny * nz * u, 1);
      } else {
        f = nz * u + 0.5 * (nx + ny);
      }
    } else {
      if (6. * nx * ny * nz * u < -cube(nz) + 3. * sqr(nz) * (nx + ny) -
                                      3. * nz * (sqr(nx) + sqr(ny)) + cube(nx) +
                                      cube(ny)) {
        f = SolveCubic(
            1., -3. * (nx + ny), 3. * (sqr(nx) + sqr(ny)),
            -(cube(nx) + cube(ny)) + 6. * nx * ny * nz * u, 1);
      } else {
        f = SolveCubic(
            2., -3. * (nx + ny + nz), 3. * (sqr(nx) + sqr(ny) + sqr(nz)),
            -(cube(nx) + cube(ny) + cube(nz)) + 6. * nx * ny * nz * u, 1);
      }
    }

    return f - 0.5 * (nx + ny + nz);
  }
  static Scal GetLineA0(Vect4 n, Scal u) {
    (void)n;
    (void)u;
    return 0; // XXX not implemented
  }
  // Plane constant by volume fraction in unit cell.
  // n : normal
  // u: volume fraction
  // Returns:
  // a: plane constant
  // Equation of reconstructed line
  // x.dot(n) = a
  template <class Vect>
  static Scal GetLineA1(const Vect& n0, Scal u) {
    Vect n = n0.abs();
    Sort(n);
    u = Clip(u, 0, 1);
    if (u <= 0.5) {
      return GetLineA0(n, u);
    }
    return -GetLineA0(n, 1 - u);
  }

  // Line constant by volume fraction in rectangular cell.
  // n : normal
  // u: volume fraction
  // h: cell size
  // Returns:
  // a: line constant
  // Equation of reconstructed line
  // x.dot(n) = a
  template <class Vect>
  static Scal GetLineA(const Vect& n, Scal u, const Vect& h) {
    return GetLineA1(n * h, u);
  }

  // GetLineVol() helper
  // assume dx > 0
  template <class Vect>
  static Scal GetLineVol0(
      const Vect& n, Scal a, const Vect& h, Scal dx, size_t d) {
    // Acceptor is a rectangular box adjacent to current cell.
    Vect hh = h; // acceptor size
    hh[d] = dx;
    // Line constant for line advected by dx
    // with origin at acceptor center
    // (e.g. shift 0 if dx=h[0], shift h[0]*0.5 if dx=0)
    Vect dc(0); // shift of center
    dc[d] = (h[d] - dx) * 0.5;
    const Scal aa = a - n.dot(dc); // new line constant
    const Scal uu = GetLineU(n, aa, hh); // volume fraction
    const Scal vv = hh.prod(); // acceptor volume
    Scal r = uu * vv; // result
    r = std::min(r, GetLineU(n, a, h) * h.prod()); // limit by fluid in cell
    return r;
  }

  // Volume surplus in downwind adjacent cell after advection in direction d.
  // (right if dx > 0, left if dx < 0)
  // n : normal
  // a: line constant
  // h: cell size
  // dx: advection displacement in x
  // d: direction 0,1,2
  // Returns:
  // volume surplus in adjacent cell
  template <class Vect>
  static Scal GetLineVol(Vect n, Scal a, const Vect& h, Scal dx, size_t d) {
    if (dx < 0.) {
      n[d] = -n[d];
      dx = -dx;
    }
    return GetLineVol0(n, a, h, dx, d);
  }

  // GetLineVolStr() helper
  // assume dx > 0
  template <class Vect>
  static Scal GetLineVolStr0(
      Vect n, Scal a, Vect h, Scal dx, Scal dxu, size_t d) {
    const Scal u = GetLineU(n, a, h); // volume fraction
    Vect sh = h; // stretched size
    sh[d] = h[d] + dx - dxu;
    const Vect sn = n / sh; // stretched normal
    const Scal sa = GetLineA(sn, u, sh); // stretched line constant
    Vect dc(0); // shift of center
    dc[d] = dxu;
    return GetLineVol0(sn, sa, sh, dx, d);
  }

  // Volume surplus in downwind adjacent cell after stretching in x
  // (right if dx > 0, left if dx < 0)
  // n : normal
  // a: line constant
  // h: cell size
  // dx: displacement in x
  // dxu: displacement in x of the other face upwind
  // d: direction 0,1,2
  // Returns:
  // volume surplus in adjacent cell
  template <class Vect>
  static Scal GetLineVolStr(
      Vect n, Scal a, Vect h, Scal dx, Scal dxu, size_t d) {
    if (dx < 0.) {
      n[d] = -n[d];
      dx = -dx;
      dxu = -dxu;
    }
    return GetLineVolStr0(n, a, h, dx, dxu, d);
  }

  // Fluid volume flux to downwind adjacent cell in x.
  // n : normal
  // h: cell size
  // q: mixture volume flux
  // dt: time step
  // d: direction 0,1,2
  // Returns:
  // fluid volume flux
  template <class Vect>
  static Scal GetLineFlux(Vect n, Scal a, Vect h, Scal q, Scal dt, size_t d) {
    const Scal s = h.prod() / h[d]; // face area
    const Scal dx = q / s * dt; // displacement
    Scal v = GetLineVol(n, a, h, dx, d);
    if (q < 0.) {
      v = -v;
    }
    return v / dt;
  }

  // Fluid volume flux to downwind adjacent cell in x after stretching.
  // n : normal
  // h: cell size
  // q: mixture volume flux
  // qu: mixture volume flux upwind face
  // dt: time step
  // d: direction 0,1,2
  // Returns:
  // fluid volume flux
  template <class Vect>
  static Scal GetLineFluxStr(
      Vect n, Scal a, Vect h, Scal q, Scal qu, Scal dt, size_t d) {
    const Scal s = h.prod() / h[d]; // face area
    const Scal dx = q / s * dt; // displacement
    const Scal dxu = qu / s * dt; // displacement on upwind face
    Scal v = GetLineVolStr(n, a, h, dx, dxu, d);
    if (q < 0.) {
      v = -v;
    }
    return v / dt;
  }

  // Returns projection of point to to plane 'n.dot(x) = a'
  // x: target point
  // n: normal
  // a: constant
  template <class Vect>
  static Vect GetPlaneProj(Vect x, Vect n, Scal a) {
    return x - n * ((n.dot(x) - a) / n.sqrnorm());
  }

  // Solves 3x3 singular system Ax=0
  // a: matrix a of rank 2, a[i] is row (equation)
  // Returns:
  // x: solution
  static Vect3 SolveSingular(const generic::Vect<Vect3, 3>& a) {
    auto p = [](size_t i) { return (i + 1) % 3; };
    auto pp = [](size_t i) { return (i + 2) % 3; };

    // d[i][j] is det of matrix without row i and column j
    generic::Vect<Vect3, 3> d;
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        d[i][j] =
            a[p(i)][p(j)] * a[pp(i)][pp(j)] - a[p(i)][pp(j)] * a[pp(i)][p(j)];
      }
    }

    // mi[j] = argmax_i(d[i][j])
    generic::Vect<size_t, 3> mi(0);
    for (size_t j = 0; j < 3; ++j) {
      for (size_t i = 1; i < 3; ++i) {
        if (std::abs(d[i][j]) > std::abs(d[mi[j]][j])) {
          mi[j] = i;
        }
      }
    }

    // j = argmax_j(d[mi[j]][j])
    size_t j = 0;
    for (size_t jj = 1; jj < 3; ++jj) {
      if (std::abs(d[mi[jj]][jj]) > std::abs(d[mi[j]][j])) {
        j = jj;
      }
    }

    size_t i = mi[j];

    Vect3 x;
    x[j] = 1.;

    // solve for other from given x[mj]
    x[p(j)] = -(a[p(i)][j] * a[pp(i)][pp(j)] - a[p(i)][pp(j)] * a[pp(i)][j]) *
              x[j] / d[i][j];
    x[pp(j)] = -(a[p(i)][p(j)] * a[pp(i)][j] - a[p(i)][j] * a[pp(i)][p(j)]) *
               x[j] / d[i][j];

    return x;
  }

  // Normal of fitting plane with equation n.dot(x) = a.
  // xx: points
  // Returns:
  // n: normal
  template <class Vect>
  static Vect GetFitN(std::vector<Vect> xx) {
    auto xc = GetCenter(xx);
    for (auto& x : xx) {
      x -= xc;
    }

    generic::Vect<Vect, 3> a;
    for (size_t i = 0; i < 3; ++i) {
      a[i] = Vect(0);
      for (auto& x : xx) {
        a[i] += x * x[i];
      }
    }

    return SolveSingular(a);
  }

  // Returns polygon cut by cell, cell center at 0, unit cell,
  // n: normal, assume 0 <= nx <= ny
  // a: line constant, assume a <= 0
  static std::vector<Vect2> GetCutPoly0(const Vect2& n, Scal a) {
    std::vector<Vect2> xx; // result

    const Scal f = a + 0.5 * n.sum();
    const Scal nx = n[0], ny = n[1];
    const Vect2 b(-0.5); // base

    // equation of plane:
    // x.dot(n) = f
    // where x=(0,0) is lower corner and x=(1,1) is upper corner
    //  y______
    //  |      |
    //  |      |
    //  |______|x
    // conditions:
    // a <= 0
    // f <= 0.5 * n.sum()
    // nx + ny >= 2 * f

    auto P = [&xx, &b](Scal x0, Scal x1) { xx.push_back(b + Vect2(x0, x1)); };

    if (f <= 0) {
      P(0, 0);
      P(0, 0);
    } else if (nx >= f) {
      P(f / nx, 0);
      P(0, f / ny);
    } else {
      P(0, f / ny);
      P(1, (f - nx) / ny);
    }

    return xx;
  }
  // Returns polygon cut by cell, cell center at 0, unit cell,
  // n: normal, assume 0 <= nx <= ny <= nz
  // a: line constant, assume a <= 0
  static std::vector<Vect3> GetCutPoly0(const Vect3& n, Scal a) {
    std::vector<Vect3> xx; // result

    Scal f = a + 0.5 * n.sum();
    Scal nx = n[0], ny = n[1], nz = n[2];
    Vect3 b(-0.5); // base

    // equation of plane:
    // x.dot(n) = f
    // where x=(0,0,0) is lower corner and x=(1,1,1) is upper corner
    //    z|                  z______
    //     |                 /|      |
    //     |_____y          / |      |
    //    /                /  |      |
    //   /x                |  |______|y
    //                     | /      /
    //                     |/______/
    //                     x
    // conditions:
    // a <= 0
    // f <= 0.5 * n.sum()
    // nx + ny + nz >= 2 * f

    auto P = [&xx, &b](Scal x0, Scal x1, Scal x2) {
      xx.push_back(b + Vect3(x0, x1, x2));
    };

    if (f <= 0) {
      P(0., 0., 0.);
      P(0., 0., 0.);
      P(0., 0., 0.);
    } else if (nx >= f) {
      P(f / nx, 0., 0.);
      P(0., f / ny, 0.);
      P(0., 0., f / nz);
    } else if (ny >= f) {
      P(1., 0., (f - nx) / nz);
      P(1., (f - nx) / ny, 0.);
      P(0., f / ny, 0.);
      P(0., 0., f / nz);
    } else if (nz >= f) {
      if (nx + ny >= f) { // ny>=f/2, ny<f, nx>=f-ny, nx>0
        P(1., 0., (f - nx) / nz);
        P(1., (f - nx) / ny, 0.);
        P((f - ny) / nx, 1., 0.);
        P(0., 1., (f - ny) / nz);
        P(0., 0., f / nz);
      } else { // nx+ny<f
        P(0., 0., f / nz);
        P(1., 0., (f - nx) / nz);
        P(1., 1., (f - nx - ny) / nz);
        P(0., 1., (f - ny) / nz);
      }
    } else {
      P(1., 0., (f - nx) / nz);
      P(1., (f - nx) / ny, 0.);
      P((f - ny) / nx, 1., 0.);
      P(0., 1., (f - ny) / nz);
      P(0., (f - nz) / ny, 1.);
      P((f - nz) / nx, 0., 1.);
    }

    return xx;
  }
  // Returns polygon cut by cell, cell center at 0, unit cell,
  // n: normal, assume 0 <= nx <= ny <= nz
  // a: line constant, assume a <= 0
  static std::vector<Vect4> GetCutPoly0(const Vect4& n, Scal a) {
    (void)n;
    (void)a;
    return {}; // XXX not implemented
  }

  // Returns polygon cut by cell, cell center at 0, unit cell
  // n0: normal
  // a: line constant
  template <class Vect>
  static std::vector<Vect> GetCutPoly1(Vect n0, Scal a) {
    if (a > 0) {
      n0 = -n0;
      a = -a;
    }
    auto n = n0.abs();
    auto r = Argsort(n);
    Vect nsorted;
    for (size_t d = 0; d < Vect::dim; ++d) {
      nsorted[d] = n[r[d]];
    }
    auto xx = GetCutPoly0(nsorted, a);
    for (auto& x : xx) {
      Vect t = x;
      for (size_t d = 0; d < Vect::dim; ++d) {
        x[r[d]] = t[d];
      }
      for (size_t d = 0; d < Vect::dim; ++d) {
        if (n0[d] < 0.) {
          x[d] *= -1.;
        }
      }
    }
    return xx;
  }

  // Returns polygon cut by cell, cell center at 0
  // n: normal
  // a: line constant
  // h: cell size
  template <class Vect>
  static std::vector<Vect> GetCutPoly2(const Vect& n, Scal a, const Vect& h) {
    auto xx = GetCutPoly1(n * h, a);
    for (auto& x : xx) {
      x *= h;
    }
    return xx;
  }

  // Returns polygon cut by cell
  // xc: cell center
  // n: normal
  // a: line constant
  // h: cell size
  template <class Vect>
  static std::vector<Vect> GetCutPoly(
      const Vect& xc, const Vect& n, Scal a, const Vect& h) {
    auto xx = GetCutPoly2(n, a, h);
    for (auto& x : xx) {
      x += xc;
    }
    return xx;
  }

  // Returns center (mean) of polygon cut by cell centered at 0
  // n: normal
  // a: line constant
  // h: cell size
  template <class Vect>
  static Vect GetCenter(const Vect& n, Scal a, const Vect& h) {
    auto xx = GetCutPoly2(n, a, h);
    return GetCenter(xx);
  }

  // Check if two points (or their projections) lie on
  // the same side from a line on plane.
  // x,xs: target points
  // x0,x1: two points defining line
  // n: plane normal
  // Returns true if:
  //   t.dot(xp-xc) * t.dot(xsp - xc) >= 0
  //   where t = n.cross(x1 - x0), xc = 0.5 * (x0 + x1)
  template <class Vect>
  static bool IsSameSide(
      const Vect& x, const Vect& xs, const Vect& x0, const Vect& x1,
      const Vect& n) {
    Vect xc = (x0 + x1) * 0.5;
    Vect t = n.cross(x1 - x0);
    return (t.dot(x - xc) >= 0.) == (t.dot(xs - xc) >= 0.);
  }

  // Check if projection of point lies inside convex polygon on plane.
  // x: target point
  // xx: polygon points
  // xc: polygon center
  // n: normal
  template <class Vect>
  static bool IsInside(
      const Vect& x, const std::vector<Vect>& xx, const Vect& xc,
      const Vect& n) {
    size_t s = xx.size();
    for (size_t i = 0; i < s; ++i) {
      if (!IsSameSide(x, xc, xx[i], xx[(i + 1) % s], n)) {
        return false;
      }
    }
    return true;
  }

  // Nearest point to convex polygon on plane.
  // x: target point
  // xx: polygon points
  // n: normal
  template <class Vect>
  static Vect GetNearest(
      const Vect& x, const std::vector<Vect>& xx, const Vect& n) {
    Vect xc = GetCenter(xx);

    if (IsInside(x, xx, xc, n)) {
      return GetPlaneProj(x, n, n.dot(xc));
    }

    Vect xn; // point with minimal distance
    Scal dn = std::numeric_limits<Scal>::max(); // minimal sqrdist
    size_t s = xx.size();
    for (size_t i = 0; i < s; ++i) {
      // nearest point to edge
      Vect xe = GetNearest(x, xx[i], xx[(i + 1) % s]);
      Scal de = xe.sqrdist(x);
      if (de < dn) {
        xn = xe;
        dn = de;
      }
    }
    return xn;
  }

  // Nearest point to polygon cut by cell.
  // x: target point
  // n: normal
  // a: line constant, relative to x=0
  // h: cell size
  template <class Vect>
  static Vect GetNearest(const Vect& x, const Vect& n, Scal a, const Vect& h) {
    auto xx = GetCutPoly2(n, a, h);
    return GetNearest(x, xx, n);
  }

  // Intersection of plane convex polygon and plane.
  // xx: points of polygon
  // sx: size of xx
  // xc: point on plane
  // n: plane normal
  // Output:
  // e: line ends
  // Returns:
  // 1: non-empty intersection
  template <class Vect>
  static bool GetInterPoly(
      const Vect* xx, size_t sx, Vect xc, Vect n, std::array<Vect, 2>& e) {
    size_t j = 0; // index in e

    for (size_t i = 0; i < sx; ++i) {
      size_t ip = (i + 1) % sx;
      Vect x0 = xx[i];
      Vect x1 = xx[ip];
      if ((n.dot(x0 - xc) > 0.) != (n.dot(x1 - xc) > 0.)) { // opposite sides
        Scal l = (xc - x0).dot(n) / (x1 - x0).dot(n);
        e[j++] = x0 + (x1 - x0) * l;

        if (j == 2) {
          return true;
        }
      }
    }
    return false;
  }

  // Intersection of plane convex polygon and plane.
  // xx: points of polygon
  // xc: point on plane
  // n: plane normal
  // Output:
  // e: line ends
  // Returns:
  // 1: non-empty intersection
  template <class Vect>
  static bool GetInterPoly(
      const std::vector<Vect>& xx, const Vect xc, const Vect n,
      std::array<Vect, 2>& e) {
    return GetInterPoly(xx.data(), xx.size(), xc, n, e);
  }

  // Intersection of plane segment and straight line.
  // e: line ends
  // x0: point on line
  // t: line direction
  // Output:
  // xi: intersection point on straight line e
  // Returns:
  // 1: xi lies between e[0] and e[1]
  static bool GetInterLine(
      const std::array<Vect2, 2>& e, const Vect2 x0, const Vect2 t, Vect2& xi) {
    Scal a = t.cross_third(x0 - e[0]) / t.cross_third(e[1] - e[0]);
    xi = e[0] + (e[1] - e[0]) * a;
    if (a >= 0. && a <= 1) {
      return true;
    }
    return false;
  }

  // Cut plane convex polygon by line.
  // Keep part positively oriented with e[1]-e[0].
  // xx: points of polygon
  // e: line ends
  // TODO: describe orientation requirements
  // Output:
  // points of cut polygon
  static std::vector<Vect2> GetCutPoly(
      const std::vector<Vect2>& xx, const std::array<Vect2, 2>& e) {
    Vect2 de = e[1] - e[0];
    Vect2 n(-de[1], de[0]); // normal, <de,n> positively oriented
    std::vector<Vect2> r; // result

    for (size_t i = 0; i < xx.size(); ++i) {
      size_t ip = (i + 1 == xx.size() ? 0 : i + 1);
      Vect2 x = xx[i];
      Vect2 xp = xx[ip];
      bool s = ((x - e[0]).dot(n) > 0.); // x inside
      bool sp = ((xp - e[0]).dot(n) > 0.); // xp inside
      if (s && sp) { // both inside
        r.push_back(xp);
      } else {
        if (s || sp) { // opposite sides
          Vect2 xi; // intersection between <s,sp> and e
          GetInterLine({x, xp}, e[0], de, xi);
          r.push_back(xi);
          if (sp) {
            r.push_back(xp);
          }
        }
      }
    }
    return r;
  }

  // Intersection of 2d convex polygon and line.
  // xx: points of polygon
  // sx: size of xx
  // xc: point on line
  // t: line tangent
  // Assume zero z-component.
  // Output:
  // aa: intersection points in form: xc + t*a
  // Returns:
  // 1: non-empty intersection
  static bool GetInter(
      const Vect2* xx, size_t sx, const Vect2& xc, const Vect2& t,
      std::array<Scal, 2>& aa) {
    size_t j = 0; // index in e

    // Intersection of lines x0,x1 and xc + t*a
    // [xc + t * a - x0, x1 - x0] = 0
    // [xc - x0, x1 - x0] + a * [t, x1 - x0] = 0
    // a = -[xc - x0, x1 - x0] / [t, x1 - x0]
    // Condition for points on opposite sides of line xc + t * a
    // [x0 - xc, t] * [x1 - xc, t] < 0

    for (size_t i = 0; i < sx; ++i) {
      size_t ip = (i + 1 == sx ? 0 : i + 1);
      Vect2 x0 = xx[i];
      Vect2 x1 = xx[ip];
      if ((t.cross_third(x0 - xc) > 0.) != (t.cross_third(x1 - xc) > 0.)) {
        Scal a = -(xc - x0).cross_third(x1 - x0) / t.cross_third(x1 - x0);
        aa[j++] = a;

        if (j == 2) {
          return true;
        }
      }
    }
    return false;
  }

  // Length of line segment.
  // xx: points of polygon
  // n: normal to line (ignored)
  // Returns:
  // a: area, positive if <xx[0]-xc,xx[1]-xc,n> is positively oriented
  static Scal GetArea(const std::vector<Vect2>& xx, Vect2 n) {
    (void)n;
    if (xx.size() != 2) {
      return 0;
    }
    return std::abs(xx[0].dist(xx[1]));
  }
  // Area of plane convex polygon.
  // xx: points of polygon
  // n: normal to plane
  // Returns:
  // a: area, positive if <xx[0]-xc,xx[1]-xc,n> is positively oriented
  static Scal GetArea(const std::vector<Vect3>& xx, Vect3 n) {
    size_t sx = xx.size();

    if (!sx) {
      return 0;
    }

    // unit normal
    n /= n.norm();
    // polygon center
    auto xc = GetCenter(xx);
    // direction
    auto xd = (xx[0] - xc);
    xd /= xd.norm();

    Scal a = 0;
    for (size_t i = 0; i < sx; ++i) {
      size_t ip = (i + 1 == sx ? 0 : i + 1);
      const Vect3 x0 = xx[i];
      const Vect3 x1 = xx[ip];
      const Vect3 s = (x1 - x0).cross(n); // surface element
      const Scal z = ((x0 + x1) * 0.5 - xc).dot(xd); // position along xd
      a += z * s.dot(xd);
    }
    return a;
  }
  static Scal GetArea(const std::vector<Vect4>& xx, Vect4 n) {
    (void)xx;
    (void)n;
    return 0; // XXX not implemented
  }
};
