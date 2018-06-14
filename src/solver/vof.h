#pragma once

#include <exception>
#include <fstream>
#include <array>
#include <memory>
#include <limits>

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
Scal SolveCubic(Scal a, Scal b, Scal c, Scal d, int k) {
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
template <class Vect>
inline Vect GetCenter(const std::vector<Vect>& xx) {
  Vect r(0); // result
  for (auto& x : xx) {
    r += x;
  }
  return r / xx.size();
}

// Nearest point to line between ends.
// x: target point
// x0,x1: line ends
template <class Scal>
inline GVect<Scal, 3> GetNearest(const GVect<Scal, 3> x,
															   const GVect<Scal, 3> x0,
															   const GVect<Scal, 3> x1) {
  using Vect = GVect<Scal, 3>;
  Vect l = x1 - x0;
  Scal k = l.dot(x - x0) / l.sqrnorm();
  Clip(k);
  return x0 + l * k;
}

// GetLineU() helper
// assuming a < 0, 0 < nx < ny
// XXX: 2d specific
template <class Scal>
inline Scal GetLineU0(Scal nx, Scal ny, Scal a) {
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

template <class Scal>
inline void Sort(GVect<Scal, 3>& v) {
  Sort(v[0], v[1], v[2]);
}

// Returns sequence of indices r such that v[r] would be sorted
template <class Scal>
inline GVect<size_t, 3> Argsort(GVect<Scal, 3> v) {
  GVect<size_t, 3> r(0, 1, 2);
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

// GetLineU() helper for unit cell
// n : normal
// a: line constant
// Returns:
// u: volume fraction
// Equation of reconstructed line 
// x.dot(n) = a
template <class Scal>
inline Scal GetLineU1(const GVect<Scal, 3>& n0, Scal a) {
  using Vect = GVect<Scal, 3>;

  Vect n = n0.abs();
  Sort(n);
  Clip(a, -0.5 * n.sum(), 0.5 * n.sum());

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
template <class Scal>
inline Scal GetLineU(const GVect<Scal, 3>& n, Scal a, 
                     const GVect<Scal, 3>& h) {
  return GetLineU1(n * h, a);
}

// GetLineA() helper
// assuming 0 < u < 0.5, 0 < nx < ny < nz
template <class Scal>
inline Scal GetLineA0(Scal nx, Scal ny, Scal nz, Scal u) {
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
inline Scal GetLineA1(const GVect<Scal, 3>& n0, Scal u) {
  using Vect = GVect<Scal, 3>;
  Vect n = n0.abs();
  Sort(n);
  Clip(u);

  if (u < 0.5) {
    return GetLineA0(n[0], n[1], n[2], u);
  } else {
    return -GetLineA0(n[0], n[1], n[2], 1. - u);
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

// GetLineVol() helper
// assume dx > 0
template <class Scal>
inline Scal GetLineVol0(const GVect<Scal, 3>& n, Scal a, 
                        const GVect<Scal, 3>& h, Scal dx, size_t d) {
  using Vect = GVect<Scal, 3>;
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
template <class Scal>
inline Scal GetLineVol(GVect<Scal, 3> n, Scal a, 
                        const GVect<Scal, 3>& h, Scal dx, size_t d) {
  if (dx < 0.) {
    n[d] = -n[d];
    dx = -dx;
  }
  return GetLineVol0(n, a, h, dx, d);
}

// GetLineVolStr() helper
// assume dx > 0
template <class Scal>
inline Scal GetLineVolStr0(const GVect<Scal, 3>& n, Scal a, 
                           const GVect<Scal, 3>& h, Scal dx, Scal dxu, 
                           size_t d) {
  using Vect = GVect<Scal, 3>;
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
template <class Scal>
inline Scal GetLineVolStr(GVect<Scal, 3> n, Scal a, 
                          const GVect<Scal, 3>& h, Scal dx, Scal dxu,
                          size_t d) {
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
template <class Scal>
inline Scal GetLineFlux(const GVect<Scal, 3>& n, Scal a, 
                        const GVect<Scal, 3>& h, Scal q, Scal dt, size_t d) {
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
template <class Scal>
inline Scal GetLineFluxStr(const GVect<Scal, 3>& n, Scal a, 
                           const GVect<Scal, 3>& h, 
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
// XXX: 2d specific
template <class Scal>
inline GVect<Scal, 3> GetLineC(const GVect<Scal, 3>& n, Scal a,
                               const GVect<Scal, 3>& h) {
  using Vect = GVect<Scal, 3>;
  std::array<Vect, 2> e = GetLineEnds(n, a, h);
  return (e[0] + e[1]) * 0.5;
}

// Closest point to cut line.
// x: target point
// n: normal
// a: line constant
// h: cell size
// XXX: 2d specific
template <class Scal>
inline GVect<Scal, 3> GetLineNearest(const GVect<Scal, 3> x,
                                     const GVect<Scal, 3>& n, Scal a,
                                     const GVect<Scal, 3>& h) {
  using Vect = GVect<Scal, 3>;
  std::array<Vect, 2> e = GetLineEnds(n, a, h);
  return GetNearest(x, e[0], e[1]);
}


// Projection to plane 'n.dot(x) = a'
// x: target point
// n: normal
// a: constant 
template <class Scal>
inline GVect<Scal, 3> GetPlaneProj(const GVect<Scal, 3> x,
                                   const GVect<Scal, 3>& n, Scal a) {
  using Vect = GVect<Scal, 3>;
  return x - n * ((n.dot(x) - a) / n.sqrnorm());
}

// Solves 3x3 singular system Ax=0
// a: matrix a of rank 2, a[i] is row (equation)
// Returns:
// x: solution
template <class Scal>
inline GVect<Scal, 3> SolveSingular(const GVect<GVect<Scal, 3>, 3>& a) {
  auto p = [](size_t i) { return (i + 1) % 3; };
  auto pp = [](size_t i) { return (i + 2) % 3; };

  using Vect = GVect<Scal, 3>;

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
template <class Scal>
inline GVect<Scal, 3> GetFitN(std::vector<GVect<Scal, 3>> xx) {
  using Vect = GVect<Scal, 3>;

  auto xc = GetCenter(xx);
  for (auto& x : xx) {
    x -= xc;
  }

  GVect<GVect<Scal, 3>, 3> a;
  for (size_t i = 0; i < 3; ++i) {
    a[i] = Vect(0);
    for (auto& x : xx) {
      a[i] += x * x[i];
    }
  }

  return SolveSingular(a);
}

// GetCutPoly() helper for unit cell, 
// assume 0 < nx < ny < nz, a < 0.
template <class Scal>
std::vector<GVect<Scal ,3>> GetCutPoly0(const GVect<Scal, 3>& n, Scal a) {
  using Vect = GVect<Scal, 3>;
  std::vector<Vect> xx; // result

  Scal f = a + 0.5 * n.sum();
  Scal nx = n[0], ny = n[1], nz = n[2];
  Vect b(-0.5); // base

  if (nx > f) {
    xx.push_back(b + Vect(f / nx, 0., 0.));
    xx.push_back(b + Vect(0., f / ny, 0.));
    xx.push_back(b + Vect(0., 0., f / nz));
  } else if (ny > f) {
    xx.push_back(b + Vect(1., 0., (f - nx) / nz));
    xx.push_back(b + Vect(1., (f - nx) / ny, 0.));
    xx.push_back(b + Vect(0., f / ny, 0.));
    xx.push_back(b + Vect(0., 0., f / nz));
  } else if (nz > f && nx + ny > f) {
    xx.push_back(b + Vect(1., 0., (f - nx) / nz));
    xx.push_back(b + Vect(1., (f - nx) / ny, 0.));
    xx.push_back(b + Vect((f - ny) / nx, 1., 0.));
    xx.push_back(b + Vect(0., 1., (f - ny) / nz));
    xx.push_back(b + Vect(0., 0., f / nz));
  } else if (nz >= f && nx + ny <= f) {
    xx.push_back(b + Vect(0., 0., f / nz));
    xx.push_back(b + Vect(1., 0., (f - nx) / nz));
    xx.push_back(b + Vect(1., 1., (a + 0.5 * (nz - nx - ny)) / nz));
    xx.push_back(b + Vect(0., 1., (f - ny) / nz));
  } else {
    xx.push_back(b + Vect(1., 0., (f - nx) / nz));
    xx.push_back(b + Vect(1., (f - nx) / ny, 0.));
    xx.push_back(b + Vect((f - ny) / nx, 1., 0.));
    xx.push_back(b + Vect(0., 1., (f - ny) / nz));
    xx.push_back(b + Vect(0., (f - nz) / ny, 1.));
    xx.push_back(b + Vect((f - nz) / nx, 0., 1.));
  }

  return xx;
}

// GetCutPoly() helper for unit cell.
template <class Scal>
std::vector<GVect<Scal ,3>> GetCutPoly1(const GVect<Scal, 3>& n0, Scal a) {
  using Vect = GVect<Scal, 3>;
  auto n = n0.abs();
  auto r = Argsort(n);
  auto xx = GetCutPoly0(Vect(n[r[0]], n[r[1]], n[r[2]]), -std::abs(a));
  for (auto& x : xx) {
    Vect t = x;
    for (size_t d = 0 ; d < Vect::dim; ++d) {
      x[r[d]] = t[d];
    }
    for (size_t d = 0 ; d < Vect::dim; ++d) {
      if (n0[d] * a > 0.) {
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
template <class Scal>
std::vector<GVect<Scal, 3>> GetCutPoly2(const GVect<Scal, 3>& n, Scal a,
                                        const GVect<Scal, 3>& h) {
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
template <class Scal>
std::vector<GVect<Scal, 3>> GetCutPoly(const GVect<Scal, 3>& xc,
                                       const GVect<Scal, 3>& n, Scal a,
                                       const GVect<Scal, 3>& h) {
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
template <class Scal>
GVect<Scal, 3> GetCenter(const GVect<Scal, 3>& n, Scal a,
                         const GVect<Scal, 3>& h) {
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
template <class Scal>
GVect<Scal, 3> GetNearestHalf(const GVect<Scal, 3>& x,
                              const GVect<Scal, 3>& x0,
                              const GVect<Scal, 3>& x1,
                              const GVect<Scal, 3>& n) {
  using Vect = GVect<Scal, 3>;

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
template <class Scal>
GVect<Scal, 3> GetNearestHalf(const GVect<Scal, 3>& x,
                              const GVect<Scal, 3>& x0,
                              const GVect<Scal, 3>& x1,
                              const GVect<Scal, 3>& xh,
                              const GVect<Scal, 3>& n) {
  using Vect = GVect<Scal, 3>;

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
template <class Scal>
bool IsSameSide(const GVect<Scal, 3>& x, const GVect<Scal, 3>& xs,
                const GVect<Scal, 3>& x0, const GVect<Scal, 3>& x1,
                const GVect<Scal, 3>& n) {
  using Vect = GVect<Scal, 3>;

  Vect xc = (x0 + x1) * 0.5;
  Vect t = n.cross(x1 - x0);
  return (t.dot(x - xc) >= 0.) == (t.dot(xs - xc) >= 0.);
}

// Check if projection of point lies inside convex polygon on plane.
// x: target point
// xx: polygon points 
// xc: polygon center
// n: normal
template <class Scal>
bool IsInside(const GVect<Scal, 3>& x,
              const std::vector<GVect<Scal, 3>>& xx,
              const GVect<Scal, 3>& xc,
              const GVect<Scal, 3>& n) {
  using Vect = GVect<Scal, 3>;

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
template <class Scal>
GVect<Scal, 3> GetNearest(const GVect<Scal, 3>& x,
                          const std::vector<GVect<Scal, 3>>& xx,
                          const GVect<Scal, 3>& n) {
  using Vect = GVect<Scal, 3>;

  Vect xc = GetCenter(xx);

  if (IsInside(x, xx, xc, n)) {
    return GetPlaneProj(x, n, n.dot(xc));
  }

  Scal dn = std::numeric_limits<Scal>::max(); // minimal sqrdist
  Vect xn; // point with minimal distance
  size_t s = xx.size();
  for (size_t i = 0; i < s; ++i) {
    // nearest point to edge
    Vect xe = GetNearest(x, xx[i], xx[(i + 1) % s]);
    Scal de = xe.sqrdist(x);
    if (de < dn) {
      dn = de;
      xn = xe;
    }
  }
  return xn;
}

// Nearest point to polygon cut by cell.
// x: target point
// n: normal
// a: line constant, relative to x=0
// h: cell size
template <class Scal>
GVect<Scal, 3> GetNearest(const GVect<Scal, 3>& x,
                          const GVect<Scal, 3>& n, Scal a,
                          const GVect<Scal, 3>& h) {
  auto xx = GetCutPoly2(n, a, h);
  return GetNearest(x, xx, n);
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
  static constexpr size_t kNp = 5; // particles in on string
  static constexpr size_t kNpp = 2; // maximum strings per cell
  FieldCell<std::array<Vect, kNp * kNpp>> fcp_; // cell list
  FieldCell<std::array<Vect, kNp * kNpp>> fcpt_; // cell list tmp
  FieldCell<std::array<Scal, kNp * kNpp>> fcpw_; // cell list weight
  FieldCell<size_t> fcps_; // cell list size

 public:
  struct Par {
    size_t dim = 3; // dimension (dim=2 assumes zero velocity in z)
    bool curvgrad = false; // compute curvature using gradient
    bool part = false; // particles
    Scal part_relax = 1.; 
    Scal part_h0 = 1.; // dist init
    Scal part_h = 1.;  // dist eq
    Scal part_kstr = 1.; // stretching
    Scal part_kattr = 1.; // attraction to reconstructed interface
    Scal part_kbend = 1.; // bending
    bool part_bendmean = true; // bending to mean angle (fit circle)
    bool part_n = false; // normal from particles
    // curvature from particles
    // if true, GetCurv returns fckp_
    bool part_k = false; // curvature from particles
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

    for (auto c : m.Cells()) {
      const Scal th = 1e-3;
      if (uc[c] > th && uc[c] < 1. - th) {
        Vect x = m.GetCenter(c) + GetCenter(fc_n_[c], fc_a_[c], h);
        Vect n = fc_n_[c];
        n /= n.norm();
        // direction in which normal has minimal component
        size_t d = n.abs().argmin(); 
        Vect xd(0); 
        xd[d] = 1.;
        // t0 orthogonal to n and d, <n,d,t0> positively oriented
        Vect t0 = n.cross(xd); 
        t0 /= t0.norm();
        // t1 orthogonal to n and t0
        Vect t1 = n.cross(t0);
        t1 /= t1.norm();
        const Scal pd = hm * par->part_h0; // distance between particles
        if (fcps_[c] == 0) { // if no particles yet
          for (int i = 0; i < kNp; ++i) {
            fcp_[c][fcps_[c]++] = x + t0 * ((Scal(i) - kNp / 2) * pd);
          }
          if (par->dim == 3) {
            for (int i = 0; i < kNp; ++i) {
              fcp_[c][fcps_[c]++] = x + t1 * ((Scal(i) - kNp / 2) * pd);
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
              << "," << c.GetRaw() << "\n";
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
  }
  void StartStep() override {
    auto sem = m.GetSem("start");
    if (sem("rotate")) {
      this->ClearIter();
      fc_u_.time_prev = fc_u_.time_curr;
      fc_u_.iter_curr = fc_u_.time_prev;
    }

    if (sem.Nested("reconst")) {
      if (this->GetTime() == 0.) {
          Reconst(fc_u_.time_curr);
      }
    }

      /*
      if (par->part) {
        if (sem("seed")) {
          SeedParticles(fc_u_.time_curr);
          if (par->dmp->Try(this->GetTime() + this->GetTimeStep(), 
                            this->GetTimeStep())) {
            DumpParticles(par->part_maxiter - 1);
          }
        }
      }
      */
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
  // Normal with height function evaluated at only two points
  void CalcNormalHeight(const FieldCell<Scal>& uc) {
    using MIdx = typename M::MIdx;
    using Dir = typename M::Dir;
    auto& bc = m.GetBlockCells();

    fc_n_.Reinit(m, Vect(0));
    fck_.Reinit(m, 0);

    for (auto c : m.SuCells()) {
      Vect bn; // best normal
      Scal bhx, bhy; // best first derivative
      Scal bk; // best curvature[k]
      // direction of plane normal
      for (Dir dn : {Dir::i, Dir::j, Dir::k}) {
        // directions of plane tangents ([d]irection [t]angents)
        Dir dtx((size_t(dn) + 1) % dim); 
        Dir dty((size_t(dn) + 2) % dim); 

        MIdx w = bc.GetMIdx(c);

        // offset in normal direction
        MIdx on = MIdx(dn);
        // offset in dtx,dty
        MIdx otx = MIdx(dtx);
        MIdx oty = MIdx(dty);
        // mesh step
        const Scal lx = m.GetCenter(c).dist(m.GetCenter(bc.GetIdx(w - otx)));
        const Scal ly = m.GetCenter(c).dist(m.GetCenter(bc.GetIdx(w - oty)));
        const Scal ln = m.GetCenter(c).dist(m.GetCenter(bc.GetIdx(w - on)));
        // evaluates height function
        // o: offset from w
        auto hh = [&](MIdx o) -> Scal {
          return 
            (uc[bc.GetIdx(w + o - on)] + 
            uc[bc.GetIdx(w + o)] + 
            uc[bc.GetIdx(w + o + on)]) * ln;
        };
        // height function
        const Scal h = hh(MIdx(0));
        const Scal hxm = hh(-otx);
        const Scal hxp = hh(otx);
        const Scal hym = hh(-oty);
        const Scal hyp = hh(oty);
        // corners: hxy
        const Scal hmm = hh(-otx - oty);
        const Scal hmp = hh(-otx + oty);
        const Scal hpm = hh(otx - oty);
        const Scal hpp = hh(otx + oty);

        // first derivative (slope)
        Scal hx = (hxp - hxm) / (2. * lx);  // centered
        Scal hy = (hyp - hym) / (2. * ly); 
        // sign: +1 if u increases in dn
        Scal sg = 
            (uc[bc.GetIdx(w + on)] - uc[bc.GetIdx(w - on)] > 0. ? 1. : -1.);
        // second derivative 
        Scal hxx = (hxp - 2. * h + hxm) / (lx * lx);
        Scal hyy = (hyp - 2. * h + hym) / (ly * ly);
        Scal hxy = ((hpp - hmp) - (hpm - hmm)) / (4. * lx * ly);
        // curvature
        //Scal k = -(hxp - 2. * h + hxm) / lx / 
        //    std::pow(1. + hx * hx, 3. / 2.); // 2d
        Scal k = (2. * hx * hy * hxy 
            -(sqr(hy) + 1.) * hxx -(sqr(hx) + 1.) * hyy) / 
            std::pow(sqr(hx) + sqr(hy) + 1., 3. / 2.);
        // outer normal
        Vect n;
        n[size_t(dtx)] = -hx;
        n[size_t(dty)] = -hy;
        n[size_t(dn)] = -sg;
        // select best with minimal slope
        if (dn == Dir::i || 
            std::abs(hx) + std::abs(hy) < std::abs(bhx) + std::abs(bhy)) {
          bn = n;
          bhx = hx;
          bhy = hy;
          bk = k;
        } 
      }
      fc_n_[c] = bn;
      Scal u = uc[c];
      const Scal th = 1e-6;
      fck_[c] = bk * (u > th && u < 1. - th ? 1. : 0.);
    }
  }
  Vect GetCellSize() const {
    Vect h; // cell size
    // XXX: specific for structured 3D mesh
    IdxCell c0(0);
    h = m.GetNode(m.GetNeighbourNode(c0, 7)) - 
        m.GetNode(m.GetNeighbourNode(c0, 0));
    assert(std::abs(h.prod() - m.GetVolume(c0)) < 1e-10);
    return h;
  }
  void Part(const FieldCell<Scal>& uc, typename M::Sem& sem) {
    if (sem("part-comma")) {
      m.Comm(&fc_n_);
      m.Comm(&fc_a_);
    }

    if (sem("part")) {
      auto& bc = m.GetBlockCells();
      auto& bn = m.GetBlockNodes();
      using MIdx = typename M::MIdx;
      Vect h = GetCellSize();
      Scal hm = h.norminf();

      SeedParticles(uc);

      const int sw = 2; // stencil halfwidth, [-sw,sw]
      const int sn = sw * 2 + 1; // stencil size
      // block offset
      GBlock<IdxCell, dim> bo(
          MIdx(-sw, -sw, par->dim == 2 ? 0 : -sw), 
          MIdx(sn, sn, par->dim == 2 ? 1 : sn)); 

      bool dm = par->dmp->Try(this->GetTime() + this->GetTimeStep(), 
                              this->GetTimeStep());
      if (dm) {
        DumpParticles(0);
      }

      // oriented angle from (xi-xim) to (xip-xi) 
      auto an = [&](int i, IdxCell c) {
        int im = i - 1;
        int ip = i + 1;
        Vect x = fcp_[c][i];
        Vect xm = fcp_[c][im];
        Vect xp = fcp_[c][ip];
        Vect dm = x - xm;
        Vect dp = xp - x;
        Scal lm = dm.norm();
        Scal lp = dp.norm();
        Scal sin = dm.cross_third(dp) / (lm * lp);
        return std::asin(sin);
      };

      // advance particles
      for (size_t it = 0; it < par->part_maxiter; ++it) {
        // compute correction
        for (auto c : m.Cells()) {
          // mean angle
          Scal anm = (an(1,c) + an(2,c) + an(3,c)) / 3.; 

          // clear force (needed for bending as applied from each corner)
          for (int i = 0; i < fcps_[c]; ++i) {
            fcpt_[c][i] = Vect(0); 
          }

          // traverse strings, assume fcps_[c] % kNp == 0
          for (int is = 0; is < fcps_[c] / kNp; ++is) {
            int i0 = is * kNp;
            int i1 = (is + 1) * kNp;
            // traverse particles, append to force
            for (int i = i0; i < i1; ++i) {
              const Vect& x = fcp_[c][i];
              Vect& t = fcpt_[c][i]; 

              // springs to adjacent particles
              for (int ii : {i - 1, i + 1}) {
                if (ii >= i0 && ii < i1) {
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
                Scal u = uc[cc];
                auto n = fc_n_[cc];
                Scal th = 1e-3;
                if (u > th && u < 1. - th) {
                  // nearest point to interface in cc
                  Vect xcc = m.GetCenter(cc);
                  Vect xn = xcc + GetNearest(x - xcc, n, fc_a_[cc], h);
                  //Vect xn = xcc + GetCenter(n, fc_a_[cc], h);
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
              if (i > i0 && i < i1 - 1) {
                int im = i - 1;
                int ip = i + 1;
                Vect xm = fcp_[c][im];
                Vect xp = fcp_[c][ip];
                Vect dm = x - xm;
                Vect dp = xp - x;
                Scal lm = dm.norm();
                Scal lp = dp.norm();

                /*
                // normals to segments, <nm,dm> positively oriented
                Vect nm = Vect(dm[1], -dm[0], 0.) / lm; 
                Vect np = Vect(dp[1], -dp[0], 0.) / lp;
                // torque
                Scal t = par->part_kbend * lm * lp * 
                  (an(i,c) - anm * (par->part_bendmean ? 1. : 0.));
                // forces
                Vect fm = nm * (t / (2. * lm));
                Vect fp = np * (t / (2. * lp));
                */
                
                Vect fm = (dm * (dm.dot(dp) / (lm * lm)) + dp) *
                    (par->part_kbend * std::sqrt(std::max(0.,
                        (lm * lp - dm.dot(dp)) / (lm * lp + dm.dot(dp)))));
                Vect fp = (dp * (-dm.dot(dp) / (lp * lp)) + dm) *
                    (par->part_kbend * std::sqrt(std::max(0.,
                        (lm * lp - dm.dot(dp)) / (lm * lp + dm.dot(dp)))));

                /*
                if (IsNan(fm) || IsNan(fp) || lm < 1e-10 || lp < 1e-10) {
                  std::cerr
                      << "c=" << c.GetRaw()
                      << " i=" << i
                      << " im=" << im
                      << " ip=" << ip
                      << " x=" << x
                      << " xm=" << xm
                      << " xp=" << xp
                      << " dm=" << dm
                      << " dp=" << dp
                      << " lm=" << lm
                      << " lp=" << lp
                      << " fm=" << fm
                      << " fp=" << fp
                      << " dm.dot(dp)=" << dm.dot(dp)
                      << std::endl;
                  std::terminate();
                }
                */

                // apply
                fcpt_[c][im] += fm;
                fcpt_[c][ip] += fp;
                fcpt_[c][i] -= (fm + fp);
              }
            }
          }
        }

        // compute error
        Scal tmax = 0.;
        Scal anmax = 0.;
        Scal anavg = 0; // aveage difference from mean angle
        size_t anavgn = 0;
        for (auto c : m.Cells()) {
          if (fcps_[c]) {
            // mean angle
            Scal anm = (an(1,c) + an(2,c) + an(3,c)) / 3.; 
            for (size_t i = 0; i < fcps_[c]; ++i) {
              fcp_[c][i] += fcpt_[c][i] * par->part_relax;
              tmax = std::max(tmax, fcpt_[c][i].norm());
              if (i > 0 && i < kNp - 1) {
                Scal e = std::abs(anm - an(i, c));
                anavg += e;
                ++anavgn;
                anmax = std::max(anmax, e);
              }
            }
          }
        }
        anavg /= anavgn;
        size_t dr = std::max<size_t>(1, 
            par->part_maxiter / par->part_report_fr);
        if (it % dr == 0 || it + 1 == par->part_maxiter) {
          std::cout 
              << "it=" << it 
              << " dxmax=" << tmax 
              << " anmax=" << anmax 
              << " anavg=" << anavg 
              << std::endl;
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
        if (fcps_[c]) { // contains particles
          // computes curvatures based on angle centered at i
          auto kk = [&](int i) -> Scal {
            int im = i - 1;
            int ip = i + 1;
            Vect x = fcp_[c][i];
            Vect xm = fcp_[c][im];
            Vect xp = fcp_[c][ip];
            Vect dm = x - xm;
            Vect dp = xp - x;
            Scal lm = dm.norm();
            Scal lp = dp.norm();
            Scal lmp = lm * lp;
            Scal k = std::sqrt(2. * (lmp - dm.dot(dp))) / lmp;
            if (dm.cross_third(dp) > 0.) {
              k = -k;
            }
            return k;
          };
          int ic = fcps_[c] / 2;
          fckp_[c] = kk(ic);
        }
      }
      m.Comm(&fckp_);

      // compute normal
      if (par->part_n) {
        for (auto c : m.Cells()) {
          if (fcps_[c]) {
            int i = 2;
            int im = i - 1;
            int ip = i + 1;
            Vect x = fcp_[c][i];
            Vect xm = fcp_[c][im];
            Vect xp = fcp_[c][ip];
            Vect dm = x - xm;
            Vect dp = xp - x;
            Scal lm = dm.norm();
            Scal lp = dp.norm();
            // normals to segments, <nm,dm> positively oriented
            Vect nm = Vect(dm[1], -dm[0], 0.) / lm; 
            Vect np = Vect(dp[1], -dp[0], 0.) / lp;
            fc_n_[c] = (nm + np) * (-0.5);
          }
        }
      }
    }
  }

  void Reconst(const FieldCell<Scal>& uc) {
    auto sem = m.GetSem("reconst");
    if (sem("height")) {
      // Reconstuction with normal from height functions
      CalcNormalHeight(uc);
      auto h = GetCellSize();
      for (auto c : m.AllCells()) {
        fc_a_[c] = GetLineA(fc_n_[c], uc[c], h);
      }
      m.Comm(&fck_);
    }

    if (par->part && par->part_n) {
      // Correction with normal from particles
      Part(uc, sem);
      if (sem("parta")) {
        auto h = GetCellSize();
        for (auto c : m.AllCells()) {
          fc_a_[c] = GetLineA(fc_n_[c], uc[c], h);
        }
      }
    }
  }
  // Print column of datafield
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
  // Print column of datafield
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
    // directions, format: {dir LE, dir EI, ...}
    std::vector<size_t> dd; 
    Scal vsc; // scaling factor for ffv, used for splitting
    if (par->dim == 3) { // 3d
      if (count_ % 3 == 0) {
        dd = {0, 1, 1, 2, 2, 0};
      } else if (count_ % 3 == 1) {
        dd = {1, 2, 2, 0, 0, 1};
      } else {
        dd = {2, 0, 0, 1, 1, 2};
      }
      vsc = 0.5;
    } else {
      if (count_ % 2 == 0) {
        dd = {0, 1};
      } else {
        dd = {1, 0};
      } 
      vsc = 1.0;
    }
    for (size_t id = 0; id < dd.size(); ++id) {
      // TODO: fluxes computed twice, consider buffer
      if (sem("adv")) {
        size_t d = dd[id]; // direction as index
        Dir md(d); // direction as Dir
        MIdx wd(md); // offset in direction d
        auto& uc = fc_u_.iter_curr;
        auto& bc = m.GetBlockCells();
        auto& bf = m.GetBlockFaces();
        auto h = GetCellSize();
        auto& ffv = *ffv_; // [f]ield [f]ace [v]olume flux
        const Scal dt = this->GetTimeStep();
        for (auto c : m.Cells()) {
          auto w = bc.GetMIdx(c);
          const Scal lc = m.GetVolume(c);
          // faces
          IdxFace fm = bf.GetIdx(w, md);
          IdxFace fp = bf.GetIdx(w + wd, md);
          // mixture volume fluxes
          const Scal vm = ffv[fm] * vsc;
          const Scal vp = ffv[fp] * vsc;
          // mixture volume cfl
          const Scal sm = vm * dt / lc;
          const Scal sp = vp * dt / lc;
          const Scal ds = sp - sm;
          // upwind cells
          IdxCell cum = m.GetNeighbourCell(fm, vm > 0. ? 0 : 1);
          IdxCell cup = m.GetNeighbourCell(fp, vp > 0. ? 0 : 1);
          if (id % 2 == 0) { // Euler Implicit
            // phase 1 volume fluxes TODO: rename mixture to q, phase 1 to smth
            Scal qm = GetLineFlux(fc_n_[cum], fc_a_[cum], h, vm, dt, d);
            Scal qp = GetLineFlux(fc_n_[cup], fc_a_[cup], h, vp, dt, d);
            // phase 1 cfl
            const Scal lm = qm * dt / lc;
            const Scal lp = qp * dt / lc;
            const Scal dl = lp - lm;
            uc[c] = (uc[c] - dl) / (1. - ds);
          } else { // Lagrange Explicit
            // upwind faces
            IdxFace fum = bf.GetIdx(vm > 0. ? w - wd : w, md);
            IdxFace fup = bf.GetIdx(vp > 0. ? w : w + wd, md);
            // upwind fluxes
            Scal vum = ffv[fum] * vsc;
            Scal vup = ffv[fup] * vsc;
            // phase 1 volume fluxes
            Scal qm = GetLineFluxStr(fc_n_[cum], fc_a_[cum], h, vm, vum, dt, d);
            Scal qp = GetLineFluxStr(fc_n_[cup], fc_a_[cup], h, vp, vup, dt, d);
            // phase 1 cfl
            const Scal lm = qm * dt / lc;
            const Scal lp = qp * dt / lc;
            const Scal dl = lp - lm;
            uc[c] = uc[c] * (1. + ds) - dl;
          }
        }

        // clip
        const Scal th = 1e-6;
        for (auto c : m.Cells()) {
          Scal& u = uc[c];
          if (u < th) {
            u = 0.;
          } else if (u > 1. - th) {
            u = 1.;
          }
        }
        m.Comm(&uc);
      }
      if (sem.Nested("reconst")) {
        Reconst(fc_u_.iter_curr);
      }
    }

    // Curvature from gradient of volume fraction
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

    if (par->part) {
      Part(fc_u_.iter_curr, sem);
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
    return par->part_k ? fckp_ : fck_;
  }
  // curvature from height function
  const FieldCell<Scal>& GetCurvH() const {
    return fck_;
  }
  // curvature from particles
  const FieldCell<Scal>& GetCurvP() const {
    return fckp_;
  }
  using P::GetField;
};

} // namespace solver
