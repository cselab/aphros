#pragma once

#include <exception>
#include <fstream>
#include <array>
#include <memory>
#include <limits>
#include <cmath>

#include "geom/vect.h"

template <class Scal, bool ba=false>
class Reconst {
 public:
  using Vect2 = GVect<Scal, 2>;
  using Vect = GVect<Scal, 3>;

  static void Clip(Scal& a, Scal l, Scal u) {
    a = std::max(l, std::min(u, a));
  }

  static void Clip(Scal& a) {
    Clip(a, 0., 1.);
  }

  static Scal GetClip(Scal a) {
    Clip(a, 0., 1.);
    return a;
  }
  static Scal GetClip(Scal a, Scal a0, Scal a1) {
    Clip(a, a0, a1);
    return a;
  }

  static Scal cube(Scal a) {
    return a * a * a;
  }

  static Scal sqr(Scal a) {
    return a * a;
  }

  static Scal SolveCubic(Scal a, Scal b, Scal c, Scal d, int k) {
    Scal p = (3. * a * c - b * b) / (3. * a * a);
    p = std::min(p, 0.);
    Scal q = (2. * cube(b) - 9. * a * b * c + 27. * a * a * d) / (27. * cube(a));
    Scal r = 3. * q * std::sqrt(-3. / p) / (2. * p);
    r = std::max(-1., std::min(1., r));
    Scal t = 2. * std::sqrt(-p / 3.) *
        std::cos(1. / 3. * std::acos(r) - 2. * M_PI * k / 3.);
    Scal x = t - b / (3. * a);
    return x;
  }

  // Center of cloud of ponts.
  static Vect GetCenter(const std::vector<Vect>& xx) {
    Vect r(0); // result
    for (auto& x : xx) {
      r += x;
    }
    return r / xx.size();
  }

  // Nearest point to line between ends.
  // x: target point
  // x0,x1: line ends
  template <class Vect>
  static Vect GetNearest(const Vect x, const Vect x0, const Vect x1) {
    Vect l = x1 - x0;
    Scal k = l.dot(x - x0) / l.sqrnorm();
    Clip(k);
    return x0 + l * k;
  }

  // GetLineU() helper
  // assuming a < 0, 0 < nx < ny
  // XXX: 2d specific
  static Scal GetLineU0(Scal nx, Scal ny, Scal a) {
    Scal a1 = 0.5 * (nx - ny);
    if (a <= a1) {
      if (nx == 0.) { // TODO revise
        return 0.;
      } else {
        return sqr(a + 0.5 * (nx + ny)) / (2. * nx * ny);
      }
    } else {
      return 0.5 + a / ny;
    }
  }

  // Sort to have a <= b <= c
  static void Sort(Scal& a, Scal& b, Scal& c) {
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

  static void Sort(Vect& v) {
    Sort(v[0], v[1], v[2]);
  }

  // Returns sequence of indices r such that v[r] would be sorted
  static GVect<size_t, 3> Argsort(Vect v) {
    GVect<size_t, 3> r(0ul, 1ul, 2ul);
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

  // Volume fraction from plane constant in unit cell with ordered normal.
  // nx,ny,nz: normal, 0 <= nx <= ny <= nz
  // a: plane constant, -0.5*n.sum() <= a <= 0
  static Scal GetLineU0(Scal nx, Scal ny, Scal nz, Scal a) {
    Scal f = a + 0.5 * (nx + ny + nz);

    if (f <= 0) {
      return 0.;
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
            std::min(1., (f - ny) / nx) * sqr(f - ny)) / (6 * ny * nz);
      }
      return (2 * f - nx - ny) / (2 * nz);
    }
    return (cube(f) - cube(f - nx) - cube(f - ny) - cube(f - nz)) /
        (6 * nx * ny * nz);
  }

  // Volume fraction from plane constant in unit cell.
  // n0: normal
  // a: plane constant
  // Returns:
  // u: volume fraction
  // Equation of reconstructed plane
  // x.dot(n) = a
  static Scal GetLineU1(const Vect& n0, Scal a) {
    Vect n = n0.abs();
    Sort(n);
    Clip(a, -0.5 * n.sum(), 0.5 * n.sum());

    if (ba) {
      return GetClip(plane_volume(n[0], n[1], n[2], a), 0, 1);
    }

    if (a < 0.) {
      return GetLineU0(n[0], n[1], n[2], a);
    } else {
      return 1. - GetLineU0(n[0], n[1], n[2], -a);
    }
  }

  // Volume fraction from line constant in rectangular cell.
  // n : normal
  // a: line constant
  // h: cell size
  // Returns:
  // u: volume fraction
  static Scal GetLineU(const Vect& n, Scal a, const Vect& h) {
    return GetLineU1(n * h, a);
  }

  // Plane constant by volume fraction in unit cell.
  // nx,ny,nz: normal, 0 <= nx <= ny <= nz
  // u: volume fraction, 0 <= u <= 0.5
  // Returns:
  // a: plane constant
  static Scal GetLineA0(Scal nx, Scal ny, Scal nz, Scal u) {
    Scal f;
    if (6. * ny * nz * u < sqr(nx)) {
      f = std::pow(6. * nx * ny * nz * u, 1. / 3.);
    } else if (6. * ny * nz * u < 3. * sqr(ny) - 3 * nx * ny + sqr(nx)) {
      f = 0.5 * nx + std::sqrt(2. * ny * nz * u - sqr(nx) / 12.);
    } else if (nz > nx + ny) {
      if (2. * nz * u < nx + ny) {
        f = SolveCubic(1., -3. * (nx + ny), 
            3. * (sqr(nx) + sqr(ny)), 
            -(cube(nx) + cube(ny)) + 6. * nx * ny * nz * u, 1);
      } else {
        f = nz * u + 0.5 * (nx + ny); 
      }
    } else {
      if (6. * nx * ny * nz * u < 
          -cube(nz) + 3. * sqr(nz) * (nx + ny) 
          -3. * nz * (sqr(nx) + sqr(ny)) + cube(nx) + cube(ny)) {
        f = SolveCubic(1., -3. * (nx + ny), 
            3. * (sqr(nx) + sqr(ny)), 
            -(cube(nx) + cube(ny)) + 6. * nx * ny * nz * u, 1);
      } else {
        f = SolveCubic(2., -3. * (nx + ny + nz),
            3. * (sqr(nx) + sqr(ny) + sqr(nz)),
            -(cube(nx) + cube(ny) + cube(nz)) + 6. * nx * ny * nz * u, 1);
      }
    }

    return f - 0.5 * (nx + ny + nz);
  }

  // GetLineA() helper
  // assuming 0 < u < 0.5, 0 < nx < ny
  // XXX: 2d specific
  static Scal GetLineA0(Scal nx, Scal ny, Scal u) {
    Scal u1 = 0.5 * nx / ny;
    if (u <= u1) {
      return -0.5 * (nx + ny) + std::sqrt(2. * nx * ny * u);
    } else {
      return ny * (u - 0.5);
    }
  }

  static double plane_alpha(Scal c, Scal nx, Scal ny, Scal nz) {
    using std::sqrt;
    using std::abs;
    using std::acos;
    using std::cos;
    using std::max;
    using std::min;
    Scal alpha;
    Scal n1x, n1y, n1z;

    n1x = abs(nx); n1y = abs(ny); n1z = abs(nz);

    Scal m1, m2, m3;
    m1 = min(n1x, n1y);
    m3 = max(n1x, n1y);
    m2 = n1z;
    if (m2 < m1) {
      Scal tmp = m1;
      m1 = m2;
      m2 = tmp;
    }
    else if (m2 > m3) {
      Scal tmp = m3;
      m3 = m2;
      m2 = tmp;
    }
    Scal m12 = m1 + m2;
    Scal pr = max(6.*m1*m2*m3, 1e-50);
    Scal V1 = m1*m1*m1/pr;
    Scal V2 = V1 + (m2 - m1)/(2.*m3), V3;
    Scal mm;
    if (m3 < m12) {
      mm = m3;
      V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;
    } else {
      mm = m12;
      V3 = mm/(2.*m3);
    }

    c = GetClip(c, 0., 1.);
    Scal ch = min(c, 1. - c);
    if (ch < V1) {
      alpha = pow (pr*ch, 1./3.);
    } else if (ch < V2) {
      alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;
    } else if (ch < V3) {
      Scal p12 = sqrt (2.*m1*m2);
      Scal q = 3.*(m12 - 2.*m3*ch)/(4.*p12);
      Scal teta = acos(GetClip(q,-1.,1.))/3.;
      Scal cs = cos(teta);
      alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;
    } else if (m12 < m3) {
      alpha = m3*ch + mm/2.;
    } else {
      Scal p = m1*(m2 + m3) + m2*m3 - 1./4., p12 = sqrt(p);
      Scal q = 3.*m1*m2*m3*(1./2. - ch)/(2.*p*p12);
      Scal teta = acos(GetClip(q,-1.,1.))/3.;
      Scal cs = cos(teta);
      alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
    }
    if (c > 1./2.) alpha = 1. - alpha;

    if (nx < 0.) {
      alpha += nx;
    }
    if (ny < 0.) {
      alpha += ny;
    }
    if (nz < 0.) {
      alpha += nz;
    }

    return alpha - 0.5 * (nx + ny + nz);
  }

  static Scal plane_volume(Scal nx, Scal ny, Scal nz, Scal alpha) {
    using std::sqrt;
    using std::abs;
    using std::max;
    using std::min;
    Scal al = alpha + (nx + ny + nz)/2. +
      max(0., -nx) + max(0., -ny) + max(0., -nz);
    if (al <= 0.) {
      return 0.;
    }
    Scal tmp = abs(nx) + abs(ny) + abs(nz);
    if (al >= tmp) {
      return 1.;
    }
    if (tmp < 1e-10) { // XXX
      return 0.;
    }
    Scal n1 = abs(nx)/tmp;
    Scal n2 = abs(ny)/tmp;
    Scal n3 = abs(nz)/tmp;
    al = max(0., min(1., al/tmp));
    Scal al0 = min(al, 1. - al);
    Scal b1 = min(n1, n2);
    Scal b3 = max(n1, n2);
    Scal b2 = n3;
    if (b2 < b1) {
      tmp = b1;
      b1 = b2;
      b2 = tmp;
    } else if (b2 > b3) {
      tmp = b3;
      b3 = b2;
      b2 = tmp;
    }
    Scal b12 = b1 + b2;
    Scal bm = min(b12, b3);
    Scal pr = max(6.*b1*b2*b3, 1e-50);
    if (al0 < b1) {
      tmp = al0*al0*al0/pr;
    } else if (al0 < b2) {
      tmp = 0.5*al0*(al0 - b1)/(b2*b3) +  b1*b1*b1/pr;
    } else if (al0 < bm) {
      tmp = (al0*al0*(3.*b12 - al0) + b1*b1*(b1 - 3.*al0) +
          b2*b2*(b2 - 3.*al0))/pr;
    } else if (b12 < b3) {
      tmp = (al0 - 0.5*bm)/b3;
    } else {
      tmp = (al0*al0*(3. - 2.*al0) + b1*b1*(b1 - 3.*al0) +
          b2*b2*(b2 - 3.*al0) + b3*b3*(b3 - 3.*al0))/pr;
    }

    Scal volume = (al <= 0.5 ? tmp : 1. - tmp);
    return GetClip(volume, 0., 1.);
  }

  // Plane constant by volume fraction in unit cell.
  // n : normal
  // u: volume fraction
  // Returns:
  // a: plane constant
  // Equation of reconstructed line 
  // x.dot(n) = a
  static Scal GetLineA1(const Vect& n0, Scal u) {
    Vect n = n0.abs();
    Sort(n);
    Clip(u);
    if (ba) {
      Scal a = plane_alpha(u, n[0], n[1], n[2]);
      return GetClip(a, -0.5 * n.sum(), 0.5 * n.sum());
    }

    if (u <= 0.5) {
      return GetLineA0(n[0], n[1], n[2], u);
    }
    return -GetLineA0(n[0], n[1], n[2], 1 - u);
  }

  // Line constant by volume fraction in rectangular cell.
  // n : normal
  // u: volume fraction
  // h: cell size
  // Returns:
  // a: line constant
  // Equation of reconstructed line 
  // x.dot(n) = a
  static Scal GetLineA(const Vect& n, Scal u, const Vect& h) {
    return GetLineA1(n * h, u);
  }

  // GetLineVol() helper
  // assume dx > 0
  static Scal GetLineVol0(const Vect& n, Scal a, 
                          const Vect& h, Scal dx, size_t d) {
    // Acceptor is a rectangular box adjacent to current cell.
    Vect hh = h; // acceptor size
    hh[d] = dx;
    // Line constant for line advected by dx 
    // with origin at acceptor center
    // (e.g. shift 0 if dx=h[0], shift h[0]*0.5 if dx=0)
    Vect dc(0); // shift of center
    dc[d] = (h[d] - dx) *  0.5; 
    Scal aa = a - n.dot(dc); // new line constant
    Scal uu = GetLineU(n, aa, hh); // volume fraction
    Scal vv = hh.prod(); // acceptor volume
    Scal r = uu * vv;  // result
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
  static Scal GetLineVol(Vect n, Scal a, const Vect& h, Scal dx, size_t d) {
    if (dx < 0.) {
      n[d] = -n[d];
      dx = -dx;
    }
    return GetLineVol0(n, a, h, dx, d);
  }

  // GetLineVolStr() helper
  // assume dx > 0
  static Scal GetLineVolStr0(const Vect& n, Scal a,
                             const Vect& h, Scal dx, Scal dxu, size_t d) {
    Scal u = GetLineU(n, a, h); // volume fraction
    Vect sh = h; // stretched size
    sh[d] = h[d] + dx - dxu;
    Vect sn = n / sh; // stretched normal
    Scal sa = GetLineA(sn, u, sh); // stretched line constant
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
  static Scal GetLineVolStr(Vect n, Scal a, const Vect& h,
                            Scal dx, Scal dxu, size_t d) {
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
  static Scal GetLineFlux(const Vect& n, Scal a, const Vect& h, 
                          Scal q, Scal dt, size_t d) {
    Scal s = h.prod() / h[d];  // face area
    Scal dx = q / s * dt; // displacement
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
  static Scal GetLineFluxStr(const Vect& n, Scal a, const Vect& h,
                             Scal q, Scal qu, Scal dt, size_t d) {
    Scal s = h.prod() / h[d];  // face area
    Scal dx = q / s * dt; // displacement
    Scal dxu = qu / s * dt; // displacement on upwind face
    Scal v = GetLineVolStr(n, a, h, dx, dxu, d);
    if (q < 0.) {
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
  // XXX: 2d specific
  static std::array<Vect, 2> GetLineEnds(
      const Vect& n, Scal a, const Vect& h) {
    // equation x.dot(n) = a;
    // (cell center is 0)
    Vect hh = h * 0.5;

    // intersection with -hh
    Vect xl((a + hh[1] * n[1]) / n[0], (a + hh[0] * n[0]) / n[1], 0.); 
    // intersection with +hh
    Vect xr((a - hh[1] * n[1]) / n[0], (a - hh[0] * n[0]) / n[1], 0.); 

    std::array<Vect, 2> e{Vect(0), Vect(0)}; // default to center
    size_t i = 0;

    if (-hh[0] <= xl[0] && xl[0] <= hh[0]) {
      e[i++] = Vect(xl[0], -hh[1], 0.);
    } 
    if (-hh[0] <= xr[0] && xr[0] <= hh[0]) {
      e[i++] = Vect(xr[0], hh[1], 0.);
    } 
    if (i < 2 && -hh[1] <= xl[1] && xl[1] <= hh[1]) {
      e[i++] = Vect(-hh[0], xl[1], 0.);
    } 
    if (i < 2 && -hh[1] <= xr[1] && xr[1] <= hh[1]) {
      e[i++] = Vect(hh[0], xr[1], 0.);
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
  // XXX: 2d specific
  static Vect GetLineC(const Vect& n, Scal a, const Vect& h) {
    std::array<Vect, 2> e = GetLineEnds(n, a, h);
    return (e[0] + e[1]) * 0.5;
  }

  // Closest point to cut line.
  // x: target point
  // n: normal
  // a: line constant
  // h: cell size
  // XXX: 2d specific
  static Vect GetLineNearest(const Vect x,
                             const Vect& n, Scal a,
                             const Vect& h) {
    std::array<Vect, 2> e = GetLineEnds(n, a, h);
    return GetNearest(x, e[0], e[1]);
  }


  // Projection to plane 'n.dot(x) = a'
  // x: target point
  // n: normal
  // a: constant 
  static Vect GetPlaneProj(const Vect x, const Vect& n, Scal a) {
    return x - n * ((n.dot(x) - a) / n.sqrnorm());
  }

  // Solves 3x3 singular system Ax=0
  // a: matrix a of rank 2, a[i] is row (equation)
  // Returns:
  // x: solution
  static Vect SolveSingular(const GVect<Vect, 3>& a) {
    auto p = [](size_t i) { return (i + 1) % 3; };
    auto pp = [](size_t i) { return (i + 2) % 3; };

    // d[i][j] is det of matrix without row i and column j
    GVect<Vect, 3> d; 
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        d[i][j] = a[p(i)][p(j)] * a[pp(i)][pp(j)] -
                  a[p(i)][pp(j)] * a[pp(i)][p(j)];
      } 
    }

    // mi[j] = argmax_i(d[i][j])
    GVect<size_t, 3> mi(0); 
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

    Vect x;
    x[j] = 1.; 

    // solve for other from given x[mj]
    x[p(j)] = -(a[p(i)][j] * a[pp(i)][pp(j)] -
               a[p(i)][pp(j)] * a[pp(i)][j]) * x[j] / d[i][j];
    x[pp(j)] = -(a[p(i)][p(j)] * a[pp(i)][j] -
               a[p(i)][j] * a[pp(i)][p(j)]) * x[j] / d[i][j];

    return x;
  }

  // Normal of fitting plane with equation n.dot(x) = a.
  // xx: points 
  // Returns:
  // n: normal
  static Vect GetFitN(std::vector<Vect> xx) {
    auto xc = GetCenter(xx);
    for (auto& x : xx) {
      x -= xc;
    }

    GVect<Vect, 3> a;
    for (size_t i = 0; i < 3; ++i) {
      a[i] = Vect(0);
      for (auto& x : xx) {
        a[i] += x * x[i];
      }
    }

    return SolveSingular(a);
  }

  // Returns polygon cut by cell, cell center at 0, unit cell,
  // n: normal, assume 0 <= nx <= ny <= nz
  // a: line constant, assume a <= 0
  static std::vector<Vect> GetCutPoly0(const Vect& n, Scal a) {
    std::vector<Vect> xx; // result

    Scal f = a + 0.5 * n.sum();
    Scal nx = n[0], ny = n[1], nz = n[2];
    Vect b(-0.5); // base

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

    auto P = [&xx,&b](Scal x0, Scal x1, Scal x2) {
      xx.push_back(b + Vect(x0, x1, x2));
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

  // Returns polygon cut by cell, cell center at 0, unit cell
  // n0: normal
  // a: line constant
  static std::vector<Vect> GetCutPoly1(Vect n0, Scal a) {
    if (a > 0) {
      n0 = -n0;
      a = -a;
    }
    auto n = n0.abs();
    auto r = Argsort(n);
    auto xx = GetCutPoly0(Vect(n[r[0]], n[r[1]], n[r[2]]), a);
    for (auto& x : xx) {
      Vect t = x;
      for (size_t d = 0 ; d < Vect::dim; ++d) {
        x[r[d]] = t[d];
      }
      for (size_t d = 0 ; d < Vect::dim; ++d) {
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
  static std::vector<Vect> GetCutPoly(const Vect& xc, const Vect& n, 
                                      Scal a, const Vect& h) {
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
  static Vect GetCenter(const Vect& n, Scal a, const Vect& h) {
    auto xx = GetCutPoly2(n, a, h);
    return GetCenter(xx);
  }

  // TODO: remove
  // Nearest point to open-ended rectangle.
  // x: target point
  // x0,x1: edge
  // n: plane normal
  // Returns:
  // xn: point on plane <xc,n> nearest to x such that
  //     t.dot(xn-xc) >= 0
  //   and projection to edge lies between x0 and x1
  //   where t = n.cross(x1-x0), xc = 0.5 * (x0 + x1)
  static Vect GetNearestHalf(const Vect& x, const Vect& x0,
                             const Vect& x1, const Vect& n) {
    Vect xc = (x0 + x1) * 0.5;
    Vect t = n.cross(x1 - x0);
    Vect l = t.cross(n);
    Vect dx = x - xc;
    Vect dxn = dx - 
        n * (n.dot(dx) / n.sqrnorm()) - 
        t * std::min(0., t.dot(dx) / t.sqrnorm()) - 
        l * std::max(0., l.dot(x - x1) / l.sqrnorm()) -
        l * std::min(0., l.dot(x - x0) / l.sqrnorm());
    return xc + dxn;
  }

  // TODO: remove
  // Nearest point to open-ended rectangle.
  // x: target point
  // x0,x1: edge of half-plane
  // xh: point on half-plane
  // n: plane normal
  // Returns:
  // xn: point on plane <xc,n> nearest to x such that
  //     t.dot(xn-xc) * t.dot(xh-xc) >= 0
  //   and projection to edge lies between x0 and x1
  //   where t = n.cross(x1-x0), xc = 0.5 * (x0 + x1)
  static Vect GetNearestHalf(const Vect& x, const Vect& x0, const Vect& x1, 
                             const Vect& xh, const Vect& n) {
    Vect xc = (x0 + x1) * 0.5;
    Vect t = n.cross(x1 - x0);
    Vect l = t.cross(n);
    Vect dx = x - xc;
    Scal sg = (t.dot(xh - xc) > 0. ? 1. : -1.);
    Vect dxn = dx - 
        n * (n.dot(dx) / n.sqrnorm()) - 
        t * (std::min(0., t.dot(dx) / t.sqrnorm() * sg) * sg) -
        l * std::max(0., l.dot(x - x1) / l.sqrnorm()) -
        l * std::min(0., l.dot(x - x0) / l.sqrnorm());
    return xc + dxn;
  }

  // Check if two points (or their projections) lie on 
  // the same side from a line on plane.
  // x,xs: target points
  // x0,x1: two points defining line
  // n: plane normal
  // Returns true if:
  //   t.dot(xp-xc) * t.dot(xsp - xc) >= 0
  //   where t = n.cross(x1 - x0), xc = 0.5 * (x0 + x1)
  static bool IsSameSide(const Vect& x, const Vect& xs,
                         const Vect& x0, const Vect& x1,
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
  static bool IsInside(const Vect& x, const std::vector<Vect>& xx,
                       const Vect& xc, const Vect& n) {
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
  static Vect GetNearest(const Vect& x, const std::vector<Vect>& xx,
                         const Vect& n) {
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
  static Vect GetNearest(const Vect& x,
                         const Vect& n, Scal a,
                         const Vect& h) {
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
  static bool GetInterPoly(const Vect* xx, size_t sx,
                           const Vect& xc, const Vect& n,
                           std::array<Vect, 2>& e) {
    size_t j = 0; // index in e

    for (size_t i = 0; i < sx; ++i) {
      size_t ip = (i + 1) % sx;
      Vect x0 = xx[i];
      Vect x1 = xx[ip];
      if ((n.dot(x0 - xc) > 0.) != (n.dot(x1 - xc) > 0.)) { // opposite sides
        Scal l = (xc - x0).dot(n)  / (x1 - x0).dot(n);
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
  static bool GetInterPoly(const std::vector<Vect>& xx, const Vect& xc,
                           const Vect& n, std::array<Vect, 2>& e) {
    return GetInterPoly(xx.data(), xx.size(), xc, n, e);
  }

  // Intersection of plane segment and straigh line.
  // e: line ends
  // x0: point on line
  // t: line direction
  // Output:
  // xi: intersection point on straight line e
  // Returns:
  // 1: xi lies between e[0] and e[1]
  static bool GetInterLine(const std::array<Vect2, 2>& e, 
                           const Vect2& x0, const Vect2& t, Vect2& xi) {
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
  static std::vector<Vect2> GetCutPoly(const std::vector<Vect2>& xx, 
                                      const std::array<Vect2, 2>& e) {
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

  // Intersection of plane and surface of fluid volume
  // xc: cell center
  // n: interface normal
  // a: line constant
  // h: cell size
  // xp: point on plane
  // np: plane normal
  // Returns:
  // nodes of polygon
  /*
  static std::vector<Vect> GetCutInter(const Vect& xc, const Vect& n, 
                                      Scal a, const Vect& h,
                                      const Vect& xp, const Vect& np) {
    // TODO
    throw std::runtime_error("not implemented");
  }
  */

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
  static bool GetInter(const Vect2* xx, size_t sx,
                       const Vect2& xc, const Vect2& t,
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

  // Area of plane convex polygon.
  // xx: points of polygon
  // n: normal to plane
  // Returns:
  // a: area, positive if <xx[0]-xc,xx[1]-xc,n> is positively oriented
  static Scal GetArea(const std::vector<Vect>& xx, Vect n) {
    size_t sx = xx.size();

    if (!sx) {
      return 0.;
    }

    // unit normal
    n /= n.norm(); 
    // polygon center
    auto xc = GetCenter(xx);
    // direction
    auto xd = (xx[0] - xc);
    xd /= xd.norm();

    Scal a = 0.;
    for (size_t i = 0; i < sx; ++i) {
      size_t ip = (i + 1 == sx ? 0 : i + 1);
      Vect x0 = xx[i];
      Vect x1 = xx[ip];
      Vect s = (x1 - x0).cross(n); // surface element
      Scal z = ((x0 + x1) * 0.5 - xc).dot(xd); // position along xd
      a += z * s.dot(xd);
    }
    return a;
  }
}; 


