#undef NDEBUG
#include <cassert>
#include <iostream>
#include <string>
#include <memory>

#include <geom/vect.h>

using Scal = double;
using Vect = GVect<Scal, 3>;

// Intersection between PLIC polygon and plane.
// xc,n,a,h: PLIC cell center, normal, plane constant, cell size
// xp,np: plane origin, normal
// Output:
// x0,x1: endpoints of intersection
// Returns true if intersect.
bool PolyInter(
    Vect xc, Vect n, Scal a, Vect h, Vect xp, Vect np, Vect& x0, Vect& x1) {
  // plane equation: np.dot(x - xp) = 0
  // PLIC plane equation: n.dot(x - xc) = a

  // plane and cell do not intersect
  if ((np * h).norm1() * 0.5 < std::abs(np.dot(xc - xp))) {
    return false;
  }
  // intersection vector
  const Vect t = n.cross(np);
  // projection of cell center on plane
  const Vect xcp = xc + np * (np.dot(xp - xc) / np.dot(np));
  // projection of cell center on PLIC
  const Vect xcl = xc + n * (a / n.dot(n));
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
  const Vect xi = xc + n * u + np * v;
  // check:
  std::cout << "n: " << n.dot(xi - xc) - a << std::endl;
  std::cout << "np: " << np.dot(xi - xp) << std::endl;
  x0 = xi;
  x1 = xi;
  return true;
}

int main() {
  Vect x0, x1;
  PolyInter(
      Vect{0., 0., 0.}, Vect{1., 0., 0.}, 0, Vect{1., 1., 1.}, Vect{0.,0.,0.}, Vect{0., 1., 0.},
      x0, x1);
  std::cout << x0 << std::endl;
}
