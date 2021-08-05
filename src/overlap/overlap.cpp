#include <float.h>
#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cmath>
using Scal = double;
#include "overlap.inc"
#include "overlap.h"

Scal GetSphereOverlap(
    const generic::Vect<Scal, 3>& x, const generic::Vect<Scal, 3>& h,
    const generic::Vect<Scal, 3>& c, Scal r) {
  Scal vv[][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
		  {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
  bool ai = true; // all inside sphere
  bool ao = true; // all outside sphere
  auto hh = h * 0.5;
  for (size_t i = 0; i < sizeof vv / sizeof *vv; ++i) {
    Scal* v = vv[i];
    v[0] = x[0] + v[0] * hh[0];
    v[1] = x[1] + v[1] * hh[1];
    v[2] = x[2] + v[2] * hh[2];
    bool in = generic::Vect<Scal, 3>(v[0], v[1], v[2]).sqrdist(c) < sqr(r);
    bool out =
	generic::Vect<Scal, 3>(v[0], v[1], v[2]).sqrdist(c) > sqr(r + h.max());
    ai = ai && in;
    ao = ao && out;
  }

  if (ai) {
    return 1.;
  }
  if (ao) {
    return 0.;
  }

  Hexahedron hex(vv[0], vv[1], vv[2], vv[3], vv[4], vv[5], vv[6], vv[7]);
  Sphere s(c[0], c[1], c[2], r);

  return overlap(s, hex) / h.prod();
}
