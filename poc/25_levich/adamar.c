#include <math.h>
#include "adamar.h"

static double sq(double x) {
  return x * x;
}

static double cube(double x) {
  return x * x * x;
}

#define DEC double mu0, mu1, a, pi, a1, a2, b1, b3, b2, U
#define DEF                                                              \
  do {                                                                   \
    mu0 = q->mu0;                                                        \
    mu1 = q->mu1;                                                        \
    a = q->a;                                                            \
    pi = q->pi;                                                          \
    a1 = pi / (a * (9 * mu1 + 6 * mu0));                                 \
    a2 = -(a * pi) / (9 * mu1 + 6 * mu0);                                \
    b1 = (sq(sq(a)) * mu1 * pi) / (9 * mu0 * mu1 + 6 * sq(mu0));         \
    b3 = (a * (2 * mu1 + 2 * mu0) * pi) / (9 * mu0 * mu1 + 6 * sq(mu0)); \
    b2 = -(sq(a) * pi) / (3 * mu0);                                      \
    U = (a * (2 * mu1 + 2 * mu0) * pi) / (9 * mu0 * mu1 + 6 * sq(mu0));  \
  } while (0)

static void f(
    struct Adamar* q, double r, double t, double* u, double* v, double* p) {
  DEC;
  DEF;
  q->U = U;
  *u = ((b2 * sq(r) + b3 * cube(r) + b1) * cos(t)) / (r * sq(r));
  *v = -((b2 * sq(r) + 2 * b3 * cube(r) - b1) * sin(t)) / (2 * r * sq(r));
  *p = (b2 * mu0 * cos(t) + pi * sq(r)) / sq(r);
}

static void g(
    struct Adamar* q, double r, double t, double* u, double* v, double* p) {
  DEC;
  DEF;
  q->U = U;
  *u = (a1 * sq(r) + a2) * cos(t);
  *v = -(2 * a1 * sq(r) + a2) * sin(t);
  *p = 10 * a1 * mu1 * r * cos(t);
}

static int c2s(
    double x, double y, double z, double* pr, double* pt, double* pp) {
  double r;
  double t;
  double p;

  r = sqrt(sq(x) + sq(y) + sq(z));
  t = (r == 0) ? 0 : acos(x / r);
  p = atan2(z, y);

  *pr = r;
  *pt = t;
  *pp = p;
  return 0;
}

static int s2c(
    double r, double t, double p, double u, double v, double* x, double* y,
    double* z) {
  *x = cos(t) * u - sin(t) * v;
  *y = cos(p) * (sin(t) * u + cos(t) * v);
  *z = sin(p) * (sin(t) * u + cos(t) * v);
  return 0;
}

int adamar_fields(
    struct Adamar* q, double x, double y, double z, double* px, double* py,
    double* pz, double* pp) {
  double r;
  double t;
  double p;
  double a;
  double u;
  double v;

  a = q->a;
  c2s(x, y, z, &r, &t, &p);
  if ((q->distance = r - a) < 0) {
    g(q, r, t, &u, &v, pp);
  } else
    f(q, r, t, &u, &v, pp);
  s2c(r, t, p, u, v, px, py, pz);
  return 0;
}
