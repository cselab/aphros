#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include <vofi.h>

const static double EPS = 1e-6;

struct {
  double x, y, z, r;
} G;

struct {
  double x, y, z, r, u, v, w;
} C;

static double sq(double x) {
  return x * x;
}

static double sphere(const double p[3]) {
  double x, y, z, r;
  x = G.x;
  y = G.y;
  z = G.z;
  r = G.r;
  return sq(p[0] - x) + sq(p[1] - y) + sq(p[2] - z) - sq(r);
}

static double circle(const double p[2]) {
  double x, y, r;
  x = G.x;
  y = G.y;
  r = G.r;
  return sq(p[0] - x) + sq(p[1] - y) - sq(r);
}

static double cylinder_distance(const double p[3]) {
  enum { X, Y, Z };
  double x, y, z, u, v, w, d;
  u = C.u;
  v = C.v;
  w = C.w;
  x = C.x - p[X];
  y = C.y - p[Y];
  z = C.z - p[Z];
  d = x * u + y * v + z * w;
  return sq(x - d * u) + sq(y - d * v) + sq(z - d * w);
}

static double cylinder(const double p[3]) {
  double r;
  r = C.r;
  return cylinder_distance(p) - sq(r);
}

double vof_3d(double x, double y, double z, double r) {
  int dim;
  double Delta, h;
  double p[3] = {-0.5, -0.5, -0.5};
  dim = 3;
  Delta = 1.0;
  G.x = x;
  G.y = y;
  G.z = z;
  G.r = r;
  h = vofi_Get_fh(sphere, NULL, Delta, dim, 0);
  return vofi_Get_cc(sphere, p, Delta, h, dim);
}

double vof_2d(double x, double y, double r) {
  int dim;
  double Delta, h;
  double p[2] = {-0.5, -0.5};
  dim = 2;
  Delta = 1.0;
  G.x = x;
  G.y = y;
  G.r = r;
  h = vofi_Get_fh(circle, NULL, Delta, dim, 0);
  return vofi_Get_cc(circle, p, Delta, h, dim);
}

double vof_cylinder(
    double x, double y, double z, double r, double u, double v, double w) {
  int dim;
  double Delta, h, d, Init;
  double p[] = {-0.5, -0.5, -0.5};
  const double* o;
  const double o0[] = {0.5, 0.5, 0.5};
  const double o1[] = {-0.5, 0.5, 0.5};
  const double o2[] = {0.5, -0.5, 0.5};

  Init = 1;
  dim = 3;
  Delta = 1.0;
  C.x = x;
  C.y = y;
  C.z = z;
  C.r = r;
  d = sq(u) + sq(v) + sq(w);
  if (r <= 0) {
    fprintf(stderr, "%s:%d: r < 0\n", __FILE__, __LINE__);
    goto err;
  }
  if (d <= 0) {
    fprintf(stderr, "%s:%d: |n| = 0\n", __FILE__, __LINE__);
    goto err;
  }
  d = sqrt(d);
  C.u = u / d;
  C.v = v / d;
  C.w = w / d;
  if (cylinder_distance(o0) > EPS)
    o = o0;
  else if (cylinder_distance(o1) > EPS)
    o = o1;
  else if (cylinder_distance(o2) > EPS)
    o = o2;
  else {
    fprintf(stderr, "%s:%d: cannot find initial guess\n", __FILE__, __LINE__);
    goto err;
  }
  h = vofi_Get_fh(cylinder, o, Delta, dim, Init);
  if (h < 0) {
    fprintf(stderr, "%s:%d: h < 0\n", __FILE__, __LINE__);
    goto err;
  }
  return vofi_Get_cc(cylinder, p, Delta, h, dim);
err:
  fprintf(
      stderr, "%s:%d: %g %g %g %g %g %g %g\n", __FILE__, __LINE__, x, y, z, r,
      u, v, w);
  return -1;
}
