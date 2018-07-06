#pragma once

#include <exception>
#include <fstream>
#include <array>
#include <memory>
#include <limits>

#include "advection.h"
#include "geom/block.h"
#include "dump/dumper.h"
#include "dump/vtk.h"

// attraction to exact sphere
#define ADHOC_ATTR 0
// normal from exact sphere
#define ADHOC_NORM 1

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
template <class Scal>
GVect<Scal, 3> GetNearest(const GVect<Scal, 3>& x,
                          const GVect<Scal, 3>& n, Scal a,
                          const GVect<Scal, 3>& h) {
  auto xx = GetCutPoly2(n, a, h);
  return GetNearest(x, xx, n);
}

// Intersection of plane convex polygon and plane.
// xx: points of polygon
// xc: point on plane
// n: plane normal
// Output:
// e: line ends
template <class Scal>
bool GetInterPoly(const std::vector<GVect<Scal, 3>>& xx,
                  const GVect<Scal, 3>& xc,
                  const GVect<Scal, 3>& n,
                  std::array<GVect<Scal, 3>, 2>& e) {
  using Vect = GVect<Scal, 3>;

  size_t j = 0; // index in e

  size_t sx = xx.size();
  for (size_t i = 0; i < xx.size(); ++i) {
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
  FieldFace<Scal> ffvu_; // flux: volume flux * field
  FieldCell<Vect> fcg_; // gradient
  FieldFace<Vect> ffg_; // gradient
  size_t count_ = 0; // number of MakeIter() calls, used for splitting
  static constexpr size_t kNp = 5; // particles in one string
  static constexpr size_t kNpp = 2; // maximum strings per cell
  FieldCell<std::array<Vect, kNp * kNpp>> fcp_; // cell list
  FieldCell<std::array<Vect, kNp * kNpp>> fcpt_; // cell list tmp
  FieldCell<std::array<Scal, kNp * kNpp>> fcpw_; // cell list weight
  FieldCell<size_t> fcps_; // cell list size

  std::vector<Vect> dpx_; // dump particles x
  std::vector<size_t> dpc_; // dump particles cell
  std::vector<std::vector<Vect>> dl_; // dump poly

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
    Scal part_intth = 1e-3; // interface detection threshold
    Scal clipth = 1e-6; // vf clipping threshold
    std::unique_ptr<Dumper> dmp; // dumper for particles
    bool dumppoly = false; // dump reconstructed interface (cut polygons)
    // XXX: adhoc
    Scal bcc_k0 = 1.;   // mul corrections to bc 
    Scal bcc_k1 = 1.;   
    Scal bcc_t0 = -1.;   // duration of phases (one negative to disable)
    Scal bcc_t1 = -1.;   
    Scal bcc_y0 = -1e10; // overwrite u=0 if y<y0 or y>y1
    Scal bcc_y1 = 1e10;  // (to remove periodic contidions)
    int part_constr = 0; // 0: no constraints
                         // 1: fixed distance, constant angle
                         // 2: fixed distance, linear angle
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
      const Scal th = par->part_intth;
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
          if (false && par->dim == 3) { // XXX adhoc no second in 3d
            for (int i = 0; i < kNp; ++i) {
              fcp_[c][fcps_[c]++] = x + t1 * ((Scal(i) - kNp / 2) * pd);
            }
          }
        }
      }
    }
  }
  // Dump cut polygons
  void DumpPoly() {
    auto sem = m.GetSem("dumppoly");
    if (sem("local")) {
      dl_.clear();
      auto h = GetCellSize();
      for (auto c : m.Cells()) {
        const Scal th = par->part_intth;
        Scal u = fc_u_.iter_curr[c];
        if (u > th && u < 1. - th) {
          dl_.push_back(GetCutPoly(m.GetCenter(c), fc_n_[c], fc_a_[c], h));
        }
      }
      using T = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<T>(&dl_));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::string st = "." + std::to_string(par->dmp->GetN());
        auto fn = "s" + st + ".vtk";
        std::cout << std::fixed << std::setprecision(8)
            << "dump" 
            << " t=" << this->GetTime() + this->GetTimeStep()
            << " to " << fn << std::endl;
        WriteVtkPoly(dl_, fn);
      }
    }
  }
  void DumpParticles(size_t it) {
    auto sem = m.GetSem("partdump");

    auto fr = par->part_dump_fr;
    size_t d = std::max<size_t>(1, par->part_maxiter / fr);
    if (fr > 1 && it % d == 0 || it + 1 == par->part_maxiter) {
      if (sem("local")) {
        dpx_.clear();
        dpc_.clear();

        // copy to arrays  
        auto& bc = m.GetBlockCells();
        for (auto c : m.Cells()) {
          for (size_t i = 0; i < fcps_[c]; ++i) {
            dpx_.push_back(fcp_[c][i]);
            auto w = bc.GetMIdx(c);
            // XXX: adhoc, hash for cell index, assume mesh size <= mn
            const size_t mn = 1000; 
            dpc_.push_back((w[2] * mn + w[1]) * mn + w[0]);
          }
        }

        // comm
        using TV = typename M::template OpCatT<Vect>;
        using TI = typename M::template OpCatT<size_t>;
        m.Reduce(std::make_shared<TV>(&dpx_));
        m.Reduce(std::make_shared<TI>(&dpc_));
      }
      if (sem("write")) {
        if (m.IsRoot()) {
          std::string st = "." + std::to_string(par->dmp->GetN());
          std::string sit = fr > 1 ? "_" + std::to_string(it) : "";
          std::string s = "partit" + st + sit + ".csv";
          std::cout << std::fixed << std::setprecision(8)
              << "dump" 
              << " t=" << this->GetTime() + this->GetTimeStep()
              << " to " << s << std::endl;
          std::ofstream o;
          o.open(s);
          o.precision(20);
          o << "x,y,z,c\n";

          for (size_t i = 0; i < dpx_.size(); ++i) {
            Vect x = dpx_[i];
            o << x[0] << "," << x[1] << "," << x[2] 
                << "," << dpc_[i] << "\n";
          }
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
  // Estimation of normal and curvature with height functions [s]
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
      std::vector<Dir> dd;
      if (par->dim == 2) {
        dd = {Dir::i, Dir::j};
      } else {
        dd = {Dir::i, Dir::j, Dir::k};
      }
      for (Dir dn : dd) {
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
        if (dn == dd[0] || 
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

      #ifdef ADHOC_NORM
      static Vect bbc;
      static Scal bbr = 0;
      static bool loaded = false;
      if (!loaded) {
        std::ifstream f("../b.dat");
        f >> bbc[0] >> bbc[1] >> bbc[2];
        f >> bbr;
        std::cout << "Loaded bbc=" << bbc << " bbr=" << bbr << std::endl;
        loaded = true;
      }
      auto q = m.GetCenter(c) - bbc;
      if (par->dim == 2) {
        q[2] = 0.;
      }
      fc_n_[c] = q / q.norm();
      #endif 
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
  // oriented angle from (x1-x0) to (x2-x1)
  static Scal GetAn(Vect x0, Vect x1, Vect x2) {
    Vect dm = x1 - x0;
    Vect dp = x2 - x1;
    Scal lm = dm.norm();
    Scal lp = dp.norm();
    Scal sin = dm.cross_third(dp) / (lm * lp);
    return std::asin(sin);
  };
  // oriented angle from (x[i]-x[i-1]) to (x[i+1]-x[i])
  static Scal GetAn(const Vect* x, size_t i) {
    return GetAn(x[i-1], x[i], x[i+1]);
  };
  // Compute force to advance particles with exact contraints on ellipse.
  // Instead of equal angles of circle, allow difference in 
  // angles -1 and +1 such that there sum equals angle 0.
  // Assume 2d positions (x,y,0).
  // x: pointer to first particle
  // sx: size of particle string
  // l: pointer to first array of line ends
  // sl: number of lines
  // par->part_kattr: relaxation factor for absolute angle
  // par->part_kbend: relaxation factor for angle between segments
  // Output:
  // f: position corrections of size sx
  void PartForce2dCE(const Vect* xx, size_t sx, 
                    const std::array<Vect, 2>* ll, size_t sl,
                    Vect* ff) {
    if (!sx) {
      return;
    }
    assert(sx == kNp && sx % 2 == 1);

    Vect h = GetCellSize();
    Scal hm = h.norminf();

    // clear force
    for (size_t i = 0; i < sx; ++i) {
      ff[i] = Vect(0);
    }

    // attracting spring to nearest point on nearest line
    for (size_t i = 0; i < sx; ++i) {
      Vect x = xx[i];

      if (sl) {
        Vect xn;  // nearest point
        Scal dn = std::numeric_limits<Scal>::max();  // sqrdist to nearest
        for (size_t j = 0; j < sl; ++j) {
          auto& e = ll[j];
          Vect xl = GetNearest(x, e[0], e[1]);
          Scal dl = xl.sqrdist(x);
          if (dl < dn) {
            xn = xl;
            dn = dl;
          }
        }

        ff[i] += (xn - x) * par->part_relax;  // scale hm
      }
    }

    // Rotates vector by pi/2
    // x: vector of plane coordinates
    auto rr = [](const Vect& x) {
      return Vect(-x[1], x[0], 0.);
    };
    // Rotates vector to angle 'a' given by unit vector
    // x: vector of plane coordinates
    // e: Vect(cos(a), sin(a), 0)
    auto re = [](const Vect& x, const Vect& e) {
      return Vect(x[0] * e[0] - x[1] * e[1], x[0] * e[1] + x[1] * e[0], 0.);
    };
    // Rotates vector to angle '-a' with 'a' given by unit vector
    // x: vector of plane coordinates
    // e: Vect(cos(a), sin(a), 0)
    auto rem = [](const Vect& x, const Vect& e) {
      return Vect(x[0] * e[0] + x[1] * e[1], -x[0] * e[1] + x[1] * e[0], 0.);
    };
    // Returns vector at angle a
    // x: vector of plane coordinates
    // e: Vect(cos(a), sin(a), 0)
    auto ra = [](Scal a) {
      return Vect(std::cos(a), std::sin(a), 0.);
    };

    // alpha: angle between x-axis and normal
    // theta: angle between segments
    // gamma: difference for angle 1
    
    // derivatives of positions by angles
    std::array<Vect, kNp> xa; // dx/dalpha
    std::array<Vect, kNp> xt; // dx/dtheta
    std::array<Vect, kNp> xg; // dx/dgamma
    // central 
    const size_t ic = (sx - 1) / 2; 
    xa[ic] = Vect(0.);
    xt[ic] = Vect(0.);
    xg[ic] = Vect(0.);
    for (size_t q = 1; q <= ic; ++q) {
      size_t i;
      Vect d;
      // forward 
      i = ic + q;
      d = xx[i] - xx[i - 1];
      xa[i] = xa[i - 1] + rr(d);
      xt[i] = xt[i - 1] + rr(d) * (q - 0.5);
      xg[i] = xg[i - 1] + rr(d) * (q - 1.);
      // backward
      i = ic - q;
      d = xx[i + 1] - xx[i];
      xa[i] = xa[i + 1] - rr(d);
      xt[i] = xt[i + 1] + rr(d) * (q - 0.5);
      xg[i] = xg[i + 1] - rr(d) * (q - 1.);
    }

    // correction of angles
    Scal da = 0.;
    Scal dt = 0.;
    Scal dg = 0.;
    for (size_t i = 0; i < sx; ++i) {
      da += ff[i].dot(xa[i]); // scale hm*hm
      dt += ff[i].dot(xt[i]);
      dg += ff[i].dot(xg[i]);
    }

    // rescale to 1
    da /= hm * hm; 
    dt /= hm * hm; 
    dg /= hm * hm; 
    
    // relaxation
    da *= par->part_kattr;
    dt *= par->part_kbend;
    dg *= par->part_kstr;

    // vector at angle da
    Vect ea = ra(da);
    // vector at angle dt/2
    Vect eth = ra(dt);
    // vector at angle dt
    Vect et = re(eth, eth);
    // vector at angle dg
    Vect eg = ra(dg);

    // segment vectors 
    std::array<Vect, kNp> dd;
    dd[ic] = Vect(0.);
    // initialize from xx
    for (size_t q = 1; q <= ic; ++q) {
      size_t i;
      // forward
      i = ic + q;
      dd[i] = xx[i] - xx[i - 1];
      // backward
      i = ic - q;
      dd[i] = xx[i + 1] - xx[i];
    }

    // apply da
    {
      for (size_t q = 1; q <= ic; ++q) {
        size_t i;
        // forward
        i = ic + q;
        dd[i] = re(dd[i], ea);
        // backward
        i = ic - q;
        dd[i] = re(dd[i], ea);
      }
    }

    // apply dt
    {
      Vect e = eth;
      for (size_t q = 1; q <= ic; ++q) {
        size_t i;
        // forward
        i = ic + q;
        dd[i] = re(dd[i], e);
        // backward
        i = ic - q;
        dd[i] = rem(dd[i], e);
        // next
        e = re(e, et);
      }
    }

    // apply dg
    {
      Vect e = eg;
      for (size_t q = 2; q <= ic; ++q) {
        size_t i;
        // forward
        i = ic + q;
        dd[i] = re(dd[i], e);
        // backward
        i = ic - q;
        dd[i] = re(dd[i], e);
        // next
        e = re(e, eg);
      }
    }

    // restore new xx from segments, store in ff
    ff[ic] = xx[ic];
    for (size_t q = 1; q <= ic; ++q) {
      size_t i;
      // forward
      i = ic + q;
      ff[i] = ff[i - 1] + dd[i];
      // backward
      i = ic - q;
      ff[i] = ff[i + 1] - dd[i];
    }

    // convert to position correction
    for (size_t i = 0; i < sx; ++i) {
      ff[i] -= xx[i];
    }
  }
  // Compute force to advance particles with exact constraints.
  // Assume 2d positions (x,y,0).
  // x: pointer to first particle
  // sx: size of particle string
  // l: pointer to first array of line ends
  // sl: number of lines
  // par->part_kattr: relaxation factor for absolute angle
  // par->part_kbend: relaxation factor for angle between segments
  // Output:
  // f: position corrections of size sx
  void PartForce2dC(const Vect* xx, size_t sx, 
                    const std::array<Vect, 2>* ll, size_t sl,
                    Vect* ff) {
    if (!sx) {
      return;
    }
    assert(sx == kNp && sx % 2 == 1);

    Vect h = GetCellSize();
    Scal hm = h.norminf();

    // clear force
    for (size_t i = 0; i < sx; ++i) {
      ff[i] = Vect(0);
    }

    // attracting spring to nearest point on nearest line
    for (size_t i = 0; i < sx; ++i) {
      Vect x = xx[i];

      if (sl) {
        Vect xn;  // nearest point
        Scal dn = std::numeric_limits<Scal>::max();  // sqrdist to nearest
        for (size_t j = 0; j < sl; ++j) {
          auto& e = ll[j];
          Vect xl = GetNearest(x, e[0], e[1]);
          Scal dl = xl.sqrdist(x);
          if (dl < dn) {
            xn = xl;
            dn = dl;
          }
        }

        if (ADHOC_ATTR) {
          auto bc = ll[sl-1][0];
          auto br = ll[sl-1][1][0];
          auto dx = x - bc;
          dx /= dx.norm();
          xn = bc + dx * br;
        }

        ff[i] += (xn - x) * par->part_relax;  // scale hm
      }
    }

    // Rotates vector by pi/2
    // x: vector of plane coordinates
    auto rr = [](const Vect& x) {
      return Vect(-x[1], x[0], 0.);
    };
    // Rotates vector to angle 'a' given by unit vector
    // x: vector of plane coordinates
    // e: Vect(cos(a), sin(a), 0)
    auto re = [](const Vect& x, const Vect& e) {
      return Vect(x[0] * e[0] - x[1] * e[1], x[0] * e[1] + x[1] * e[0], 0.);
    };
    // Rotates vector to angle '-a' with 'a' given by unit vector
    // x: vector of plane coordinates
    // e: Vect(cos(a), sin(a), 0)
    auto rem = [](const Vect& x, const Vect& e) {
      return Vect(x[0] * e[0] + x[1] * e[1], -x[0] * e[1] + x[1] * e[0], 0.);
    };
    // Returns vector at angle a
    // x: vector of plane coordinates
    // e: Vect(cos(a), sin(a), 0)
    auto ra = [](Scal a) {
      return Vect(std::cos(a), std::sin(a), 0.);
    };

    // alpha: angle between x-axis and normal
    // theta: angle between segments
    
    // derivatives of positions by angles
    std::array<Vect, kNp> xa; // dx/dalpha
    std::array<Vect, kNp> xt; // dx/dtheta
    // central 
    const size_t ic = (sx - 1) / 2; 
    xa[ic] = Vect(0.);
    xt[ic] = Vect(0.);
    for (size_t q = 1; q <= ic; ++q) {
      size_t i;
      Vect d;
      // forward 
      i = ic + q;
      d = xx[i] - xx[i - 1];
      xa[i] = xa[i - 1] + rr(d);
      xt[i] = xt[i - 1] + rr(d) * (q - 0.5);
      // backward
      i = ic - q;
      d = xx[i + 1] - xx[i];
      xa[i] = xa[i + 1] - rr(d);
      xt[i] = xt[i + 1] + rr(d) * (q - 0.5);
    }

    // correction of angles
    Scal da = 0.;
    Scal dt = 0.;
    for (size_t i = 0; i < sx; ++i) {
      da += ff[i].dot(xa[i]); // scale hm*hm
      dt += ff[i].dot(xt[i]);
    }

    // rescale to 1
    da /= hm * hm; 
    dt /= hm * hm; 
    
    // relaxation
    da *= par->part_kattr;
    dt *= par->part_kbend;

    // vector at angle da
    Vect ea = ra(da);
    // vector at angle dt/2
    Vect eth = ra(dt);
    // vector at angle dt
    Vect et = re(eth, eth);

    // segment vectors 
    std::array<Vect, kNp> dd;
    dd[ic] = Vect(0.);
    // initialize from xx
    for (size_t q = 1; q <= ic; ++q) {
      size_t i;
      // forward
      i = ic + q;
      dd[i] = xx[i] - xx[i - 1];
      // backward
      i = ic - q;
      dd[i] = xx[i + 1] - xx[i];
    }

    // apply da
    {
      for (size_t q = 1; q <= ic; ++q) {
        size_t i;
        // forward
        i = ic + q;
        dd[i] = re(dd[i], ea);
        // backward
        i = ic - q;
        dd[i] = re(dd[i], ea);
      }
    }

    // apply dt
    {
      Vect e = eth;
      for (size_t q = 1; q <= ic; ++q) {
        size_t i;
        // forward
        i = ic + q;
        dd[i] = re(dd[i], e);
        // backward
        i = ic - q;
        dd[i] = rem(dd[i], e);
        // next
        e = re(e, et);
      }
    }

    // displacement of center
    Vect dx(0);
    for (size_t i = 0; i < sx; ++i) {
      dx += ff[i];
    }
    dx *= par->part_kstr;

    // restore new xx from segments, store in ff
    ff[ic] = xx[ic];
    for (size_t q = 1; q <= ic; ++q) {
      size_t i;
      // forward
      i = ic + q;
      ff[i] = ff[i - 1] + dd[i];
      // backward
      i = ic - q;
      ff[i] = ff[i + 1] - dd[i];
    }

    // apply dx
    for (size_t i = 0; i < sx; ++i) {
      ff[i] += dx;
    }

    // convert to position correction
    for (size_t i = 0; i < sx; ++i) {
      ff[i] -= xx[i];
    }
  }
  // Compute force to advance particles.
  // Assume 2d positions (x,y,0).
  // x: pointer to first particle
  // sx: size of particle string
  // l: pointer to first array of line ends
  // nl: number of lines
  // Output:
  // f: position corrections of size sx
  void PartForce2d(const Vect* xx, size_t sx, 
                   const std::array<Vect, 2>* ll, size_t sl,
                   Vect* ff) {
    if (!sx) {
      return;
    }

    Vect h = GetCellSize();
    Scal hm = h.norminf();

    // mean angle
    Scal anm = (GetAn(xx,1) + GetAn(xx,2) + GetAn(xx,3)) / 3.; 

    // clear force
    for (size_t i = 0; i < sx; ++i) {
      ff[i] = Vect(0);
    }

    // stretching springs
    for (size_t i = 0; i < sx - 1; ++i) {
      Vect dx = xx[i + 1] - xx[i];
      Scal k = par->part_kstr;
      Scal d0 = par->part_h * hm; // equilibrium length
      Vect f = dx * (k * (1. - d0 / dx.norm()));
      ff[i] += f;
      ff[i + 1] -= f;
    }

    // attracting spring to nearest point on nearest line
    for (size_t i = 0; i < sx; ++i) {
      Vect x = xx[i];

      if (sl) {
        Vect xn;  // nearest point
        Scal dn = std::numeric_limits<Scal>::max();  // sqrdist to nearest
        for (size_t j = 0; j < sl; ++j) {
          auto& e = ll[j];
          Vect xl = GetNearest(x, e[0], e[1]);
          Scal dl = xl.sqrdist(x);
          if (dl < dn) {
            xn = xl;
            dn = dl;
          }
        }

        ff[i] += (xn - x) * par->part_kattr;
      }
    }

    // bending
    for (size_t i = 1; i < sx - 1; ++i) {
      size_t im = i - 1;
      size_t ip = i + 1;
      Vect x = xx[i];
      Vect xm = xx[im];
      Vect xp = xx[ip];
      Vect dm = x - xm;
      Vect dp = xp - x;
      Scal lm = dm.norm();
      Scal lp = dp.norm();

      // normals to segments, <nm,dm> positively oriented
      Vect nm = Vect(dm[1], -dm[0], 0.) / lm; 
      Vect np = Vect(dp[1], -dp[0], 0.) / lp;
      // torque
      Scal t = par->part_kbend * lm * lp * 
          (GetAn(xx,i) - anm * (par->part_bendmean ? 1. : 0.));
      // forces
      Vect fm = nm * (t / (2. * lm));
      Vect fp = np * (t / (2. * lp));
      // apply
      ff[im] += fm;
      ff[ip] += fp;
      ff[i] -= (fm + fp);
    }

    // scale by relaxaion factor
    for (size_t i = 0; i < sx; ++i) {
      ff[i] *= par->part_relax;
    }
  }
  // Curvature for plane string
  Scal PartK(Vect* xx, size_t sx) {
    size_t i = (sx - 1) / 2;
    Vect x = xx[i];
    Vect xm = xx[i - 1];
    Vect xp = xx[i + 1];
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
  }
  void Part(const FieldCell<Scal>& uc, typename M::Sem& sem) {
    if (sem("part-comma")) {
      m.Comm(&fc_a_);
      m.Comm(&fc_n_);
    }

    bool dm = par->dmp->Try(this->GetTime() + this->GetTimeStep(), 
                            this->GetTimeStep());

    if (sem("part-seed")) {
      SeedParticles(uc);
    }

    if (sem.Nested("part-dump0")) {
      if (dm) {
        DumpParticles(0);
      }
    }

    if (sem("part-advance")) {
      // particle strings
      std::vector<IdxCell> sc; // cell
      std::vector<Vect> xx; // particle positions
      std::vector<size_t> sx; // xx index, plus element sx.size() 
      std::vector<std::array<Vect, 2>> ll; // interface lines
      std::vector<size_t> sl; // sl index, plus element sl.size() 
      std::vector<Vect> smx; // local unit x
      std::vector<Vect> smy; // local unit y
      std::vector<Vect> smc; // local center
      std::vector<Scal> sk; // curvature

      // Extract interface, project particles
      for (auto c : m.Cells()) {
        // XXX: assume fcps_[c] % kNp == 0
        for (int is = 0; is < fcps_[c] / kNp; ++is) {
          int i0 = is * kNp;
          int i1 = (is + 1) * kNp;

          sc.push_back(c);

          // Plane coordinates
          // center
          Vect rc = fcp_[c][(i0 + i1 + 1) / 2];
          // tangent, assume straight line
          Vect rt = (fcp_[c][i0 + 1] - fcp_[c][i0]);
          rt /= rt.norm();
          // interface normal, assume orthogonal to rt
          Vect n = fc_n_[c];
          n /= n.norm();
          // string plane normal
          Vect rn = rt.cross(n);
          rn /= rn.norm();

          Vect mx = rt;
          Vect my = n;
          Vect mc = rc;

          smx.push_back(mx);
          smy.push_back(my);
          smc.push_back(mc);

          // Transform to plane coordinates.
          // x: space coordinates
          // Returns:
          // Vect(xl, yl, 0): plane coordinates
          auto pr = [&](Vect x) -> Vect {
            Vect q(0);
            q[0] = (x - mc).dot(mx);
            q[1] = (x - mc).dot(my);
            return q;
          };

          // Copy projected particles
          sx.push_back(xx.size());
          for (int i = i0; i < i1; ++i) {
            xx.push_back(pr(fcp_[c][i]));
          }

          // Extract interface lines
          sl.push_back(ll.size());
          auto& bc = m.GetBlockCells();
          using MIdx = typename M::MIdx;
          Vect h = GetCellSize();

          const int sw = 2; // stencil halfwidth, [-sw,sw]
          const int sn = sw * 2 + 1; // stencil size

          // block of offsets
          GBlock<IdxCell, dim> bo(
              MIdx(-sw, -sw, par->dim == 2 ? 0 : -sw), 
              MIdx(sn, sn, par->dim == 2 ? 1 : sn)); 

          #ifdef ADHOC_ATTR
          static Vect bbc;
          static Scal bbr = 0;
          static bool loaded = false;
          if (!loaded) {
            std::ifstream f("../b.dat");
            f >> bbc[0] >> bbc[1] >> bbc[2];
            f >> bbr;
            std::cout << "Loaded bbc=" << bbc << " bbr=" << bbr << std::endl;
            loaded = true;
          }
          #endif

          auto w = bc.GetMIdx(c);
          for (auto wo : bo) {
            auto cc = bc.GetIdx(w + wo);
            Scal u = uc[cc];
            Scal th = par->part_intth;
            if (u > th && u < 1. - th) {
              auto xcc = m.GetCenter(cc);
              auto xx = GetCutPoly(xcc, fc_n_[cc], fc_a_[cc], h);
              std::array<Vect, 2> e;
              if (GetInterPoly(xx, rc, rn, e)) {
                // projected line ends
                ll.push_back({pr(e[0]), pr(e[1])});
              }
              #ifdef ADHOC_ATTR
              ll.push_back({pr(bbc), Vect(bbr,bbr,bbr)}); 
              #endif
            }
          }
        }
      }
      sx.push_back(xx.size());
      sl.push_back(ll.size());

      assert(sx.size() == sc.size() + 1);
      assert(sl.size() == sc.size() + 1);
      assert(smx.size() == sc.size());
      assert(smy.size() == sc.size());
      assert(smc.size() == sc.size());

      sk.resize(sc.size(), 0.);

      std::vector<Vect> ff(xx.size()); // force

      // advance particles
      for (size_t it = 0; it < par->part_maxiter; ++it) {
        for (size_t i = 0; i < sc.size(); ++i) {
          if (par->part_constr == 1) {
            PartForce2dC(&(xx[sx[i]]), sx[i+1] - sx[i],
                        &(ll[sl[i]]), sl[i+1] - sl[i],
                        &(ff[sx[i]]));
          } else if (par->part_constr == 2) {
            PartForce2dCE(&(xx[sx[i]]), sx[i+1] - sx[i],
                        &(ll[sl[i]]), sl[i+1] - sl[i],
                        &(ff[sx[i]]));
          } else {
            PartForce2d(&(xx[sx[i]]), sx[i+1] - sx[i],
                        &(ll[sl[i]]), sl[i+1] - sl[i],
                        &(ff[sx[i]]));
            // freeze central particle
            Vect f = ff[(sx[i] + sx[i+1] - 1) / 2];
            for (size_t j = sx[i]; j < sx[i+1]; ++j) {
              ff[j] -= f;
            }
          }
        }

        // report error
        size_t dr = std::max<size_t>(1, 
            par->part_maxiter / par->part_report_fr);
        if (m.IsRoot() && (it % dr == 0 || it + 1 == par->part_maxiter)) {
          Scal tmax = 0.;
          Scal anmax = 0.;
          Scal anavg = 0; // average difference from mean angle
          size_t anavgn = 0;
          for (size_t i = 0; i < sc.size(); ++i) {
            const Vect* px = &(xx[sx[i]]);
            // maximum force
            for (size_t j = sx[i]; j < sx[i+1]; ++j) {
              tmax = std::max(tmax, ff[j].norm());
            }
            Scal anm = 0.;
            size_t anmn = 0;
            // mean angle in string
            for (size_t j = sx[i] + 1; j + 1 < sx[i+1]; ++j) {
              anm += GetAn(xx.data(), j);
              ++anmn;
            }
            anm /= anmn;
            // error in angle
            for (size_t j = sx[i] + 1; j + 1 < sx[i+1]; ++j) {
              Scal e = std::abs(anm - GetAn(xx.data(), j));
              anavg += e;
              ++anavgn;
              anmax = std::max(anmax, e);
            }
          }
          anavg /= anavgn;
          std::cout << std::setprecision(10)
              << "it=" << it 
              << " dxmax=" << tmax 
              << " anmax=" << anmax 
              << " anavg=" << anavg 
              << std::endl;
        }

        // advance
        for (size_t i = 0; i < xx.size(); ++i) {
          xx[i] += ff[i];
        }

        // copy back to field
        if (dm || it + 1 == par->part_maxiter) {
          fcps_.Reinit(m, 0);

          for (size_t i = 0; i < sc.size(); ++i) {
            IdxCell c = sc[i];
            Vect mx = smx[i];
            Vect my = smy[i];
            Vect mc = smc[i];
            for (size_t j = sx[i]; j < sx[i+1]; ++j) {
              auto q = xx[j];
              fcp_[c][fcps_[c]++] = mc + mx * q[0] + my * q[1];
            }
          }
        }

        if (dm) {
          //// XXX: disable iter dump to avoid suspender loop
          //DumpParticles(it + 1); 
        }
      }

      // compute curvature on strings
      for (size_t i = 0; i < sc.size(); ++i) {
        sk[i] = PartK(&(xx[sx[i]]), sx[i + 1] - sx[i]);
      }

      // compute curvature in cells
      fckp_.Reinit(m, 0.);
      {
        size_t i = 0;
        while (i < sc.size()) {
          IdxCell c = sc[i];
          Scal k = sk[i];
          ++i;
          size_t nk = 1;
          // average over all strings in c
          while (i < sc.size() && sc[i] == c) {
            k += sk[i];
            ++i;
            ++nk;
          }
          if (par->dim == 3) {
            k *= 2.;
          }
          fckp_[c] = k / nk;
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

    if (sem.Nested("part-dumplast")) {
      if (dm) {
        DumpParticles(par->part_maxiter - 1);
      }
    }
  }

  void Reconst(const FieldCell<Scal>& uc) {
    
    // XXX: adhoc
    // overwrite u=0 if y<y0 or y>y0
    {
      Scal y0 = par->bcc_y0;
      Scal y1 = par->bcc_y1;
      auto& uu = const_cast<FieldCell<Scal>&>(uc);
      for (auto c : m.AllCells()) {
        auto x = m.GetCenter(c);
        if (x[1] < y0 || x[1] > y1) {
          uu[c] = 0.;
        }
      }
    }

    auto sem = m.GetSem("reconst");
    if (sem("height")) {
      // Compute normal and curvature [s]
      CalcNormalHeight(uc);
      auto h = GetCellSize();
      // Reconstruct interface [s]
      for (auto c : m.SuCells()) {
        fc_a_[c] = GetLineA(fc_n_[c], uc[c], h);
      }
    }

    if (par->part && par->part_n) {
      Part(uc, sem);
      // Correction with normal from particles
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

        ffvu_.Reinit(m);

        for (auto f : m.Faces()) {
          auto p = bf.GetMIdxDir(f);
          MIdx wf = p.first;
          Dir df = p.second;

          if (df != md) {
            continue;
          }

          // mixture flux
          const Scal v = ffv[f] * vsc;
          // upwind cell
          IdxCell cu = m.GetNeighbourCell(f, v > 0. ? 0 : 1);
          if (id % 2 == 0) { // Euler Implicit
            // phase 2 flux
            ffvu_[f] = GetLineFlux(fc_n_[cu], fc_a_[cu], h, v, dt, d);
          } else { // Lagrange Explicit
            // upwind face
            IdxFace fu = bf.GetIdx(v > 0. ? wf - wd : wf + wd, md);
            // upwind mixture flux
            Scal vu = ffv[fu] * vsc;
            // phase 2 flux
            ffvu_[f] = GetLineFluxStr(fc_n_[cu], fc_a_[cu], h, v, vu, dt, d);
          }
        }


        ffu_.Reinit(m);
        // interpolate field value to boundaries
        InterpolateB(uc, mfc_, ffu_, m);

        // override upwind flux
        for (const auto& it : mfc_) {
          IdxFace f = it.GetIdx();
          CondFace* cb = it.GetValue().get(); 
          Scal v = ffv[f];
          if ((cb->GetNci() == 0) != (v > 0.)) {
            ffvu_[f] = v * ffu_[f];
            // XXX: adhoc
            // Alternating mul correction of flux
            // (done for bubble detachment)
            if (par->bcc_t0 > 0. && par->bcc_t1 > 0.) {
              Scal k0 = par->bcc_k0;
              Scal k1 = par->bcc_k1;
              Scal t0 = par->bcc_t0;
              Scal t1 = par->bcc_t1;
              Scal t = this->GetTime();
              Scal ts = t0 + t1;
              Scal ph = t / ts; 
              ph = ph - int(ph);
              ph *= ts;
              ffvu_[f] *= (ph < t0 ? k0 : k1);
            }
          }
        }

        for (auto c : m.Cells()) {
          auto w = bc.GetMIdx(c);
          const Scal lc = m.GetVolume(c);
          // faces
          IdxFace fm = bf.GetIdx(w, md);
          IdxFace fp = bf.GetIdx(w + wd, md);
          // mixture fluxes
          const Scal vm = ffv[fm] * vsc;
          const Scal vp = ffv[fp] * vsc;
          // mixture cfl
          const Scal sm = vm * dt / lc;
          const Scal sp = vp * dt / lc;
          const Scal ds = sp - sm;
          // phase 2 fluxes
          Scal qm = ffvu_[fm];
          Scal qp = ffvu_[fp];
          // phase 2 cfl
          const Scal lm = qm * dt / lc;
          const Scal lp = qp * dt / lc;
          const Scal dl = lp - lm;
          if (id % 2 == 0) { // Euler Implicit
            uc[c] = (uc[c] - dl) / (1. - ds);
          } else { // Lagrange Explicit
            uc[c] = uc[c] * (1. + ds) - dl;
          }
        }

        // clip
        const Scal th = par->clipth;
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
    }

    if (sem("curvcomm")) {
      m.Comm(&fck_);
    }

    if (par->dumppoly) {
      bool dm = par->dmp->Try(this->GetTime() + this->GetTimeStep(), 
                              this->GetTimeStep());
      if (dm && sem("dumppoly")) {
        DumpPoly();
      }
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
