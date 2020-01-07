#undef NDEBUG
#include <cassert>
#include <iostream>
#include <memory>
#include <string>

#include <geom/vect.h>

using Scal = double;
using Vect = GVect<Scal, 3>;

// Clips parameterized line segment.
// parametrization of segment:
// x(s) = xc + t * s
// s0,s1: endpoints
// xmin,xmax: clipping limits, xmin < xmax
// Returns true clipped segment is non-empty and updates s0 and s1.
bool ClipSegment(Scal xc, Scal t, Scal xmin, Scal xmax, Scal& s0, Scal& s1) {
  if (s0 > s1) {
    std::swap(s0, s1);
  }
  const Scal x0 = xc + t * s0;
  const Scal x1 = xc + t * s1;
  if (t >= 0) { // [x0, x1]
    if (x1 < xmin || x0 > xmax) {
      return false;
    }
    if (x0 < xmin) {
      s0 = (xmin - xc) / t;
    }
    if (x1 > xmax) {
      s1 = (xmax - xc) / t;
    }
    return true;
  }
  // t < 0, [x1, x0]
  if (x0 < xmin || x1 > xmax) {
    return false;
  }
  if (x1 < xmin) {
    s1 = (xmin - xc) / t;
  }
  if (x0 > xmax) {
    s0 = (xmax - xc) / t;
  }
  return true;
}

// Intersection between PLIC polygon and plane.
// xc,n,a,h: PLIC cell center, normal, plane constant, cell size
// xp,np: plane origin, normal
// Output:
// x0,x1: endpoints of intersection
// Returns true if intersect.
bool PolyInter(
    Vect xc, Vect n, Scal a, Vect h, Vect xp, Vect np, Vect& x0, Vect& x1) {
  const Vect hh = h * 0.5;
  // plane equation: np.dot(x - xp) = 0
  // PLIC plane equation: n.dot(x - xc) = a

  // plane and cell do not intersect
  if ((np * hh).norm1() < std::abs(np.dot(xc - xp))) {
    return false;
  }
  // intersection point witn parameters u,v:
  //   x = xc + n * u + np * v
  // intersection conditions:
  //   n.dot(x - xc) = a
  //   np.dot(x - xc) = np.dot(xp - xc)
  // intersection conditions for parameters u,v:
  //   n.dot(n) * u + n.dot(np) * v = a
  //   np.dot(n) * u + np.dot(np) * v = np.dot(xp - xc)
  // matrix:
  const Scal m11 = n.dot(n);
  const Scal m12 = n.dot(np);
  const Scal m22 = np.dot(np);
  // determinant:
  const Scal det = m11 * m22 - sqr(m12);
  if (det == 0) {
    return false;
  }
  // rhs:
  const Scal r1 = a;
  const Scal r2 = np.dot(xp - xc);
  // solution:
  const Scal u = (r1 * m22 - r2 * m12) / det;
  const Scal v = (r2 * m11 - r1 * m12) / det;
  // point on intersection line:
  const Vect xl = xc + n * u + np * v;
  const int dim = 3;
  // intersection vector
  const Vect t = n.cross(np);
  // line parametrization:
  //   x = xl + t * s
  Scal ss[2];
  const int im = t.abs().argmax();
  ss[0] = (xc[im] - hh[im] - xl[im]) / t[im];
  ss[1] = (xc[im] + hh[im] - xl[im]) / t[im];
  // clip line in other directions
  for (int i = 0; i < dim; ++i) {
    if (i == im) continue;
    if (!ClipSegment(xl[i], t[i], xc[i] - hh[i], xc[i] + hh[i], ss[0], ss[1])) {
      return false;
    }
  }
  x0 = xl + t * ss[0];
  x1 = xl + t * ss[1];
  return true;
}

int main() {
  Vect x0, x1;
  bool q = PolyInter(
      Vect{0., 0., 0.}, Vect{1., 1., 1.}, 0.3, Vect{1., 1., 1.}, Vect{0., 0., 0.},
      Vect{0., 1., 0.}, x0, x1);
  std::cout << x0 << std::endl;
  std::cout << x1 << std::endl;
  std::cout << q << std::endl;
}
