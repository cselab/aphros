#include <float.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <array>
#include <bitset>
#include <cassert>
#include <numeric>
#include <type_traits>

#include "double_prec.inc"
#include "util.inc"
#include "overlap.inc"

extern "C" {
#include "overlap.h"
}

double overlap_3d(double x, double y, double z, double r) {
  double x0, x1, y0, y1, z0, z1;

  x0 = y0 = z0 = -0.5;
  x1 = y1 = z1 = 0.5;

  vector_t v0{x0, y0, z0};
  vector_t v1{x1, y0, z0};
  vector_t v2{x1, y1, z0};
  vector_t v3{x0, y1, z0};
  vector_t v4{x0, y0, z1};
  vector_t v5{x1, y0, z1};
  vector_t v6{x1, y1, z1};
  vector_t v7{x0, y1, z1};
  Hexahedron hex{v0, v1, v2, v3, v4, v5, v6, v7};
  Sphere s{x, y, z, r};
  return overlap(s, hex);
}

double overlap_2d(double x, double y, double r) {
  double eps, x0, x1, y0, y1, z0, z1;

  eps = 0.001;
  x0 = -0.5;
  x1 = 0.5;
  y0 = -0.5;
  y1 = 0.5;
  z0 = -eps / 2;
  z1 = eps / 2;

  vector_t v0{x0, y0, z0};
  vector_t v1{x1, y0, z0};
  vector_t v2{x1, y1, z0};
  vector_t v3{x0, y1, z0};
  vector_t v4{x0, y0, z1};
  vector_t v5{x1, y0, z1};
  vector_t v6{x1, y1, z1};
  vector_t v7{x0, y1, z1};
  Hexahedron hex{v0, v1, v2, v3, v4, v5, v6, v7};
  Sphere s{x, y, 0, r};
  return overlap(s, hex) / eps;
}
