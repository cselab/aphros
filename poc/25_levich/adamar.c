#include <math.h>

struct Adamr {
  double mu0;
  double mu1;
  double a;
  double pi;
};

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

void f(struct Adamr* q) {
  DEC;
  double r;
  double t;

  double u;
  double v;
  double p;

  DEF;

  u = ((b2 * sq(r) + b3 * cube(r) + b1) * cos(t)) / (r * sq(r));
  v = -((b2 * sq(r) + 2 * b3 * cube(r) - b1) * sin(t)) / (2 * r * sq(r));
  p = (b2 * mu0 * cos(t) + pi * sq(r)) / sq(r);
}

void g(struct Adamr* q) {
  DEC;
  double r;
  double t;

  double u;
  double v;
  double p;

  DEF;
  u = (a1 * sq(r) + a2) * cos(t);
  v = -(2 * a1 * sq(r) + a2) * sin(t);
  p = 10 * a1 * mu1 * r * cos(t);
}
