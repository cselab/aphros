#include "fractions.h"
#include "curvature.h"

#include <unistd.h>


#include "learn.h"

#define SA (SW*2+1)
#define GETI(x,y) ((y)*SA+(x))
#define GET(x,y) (u[GETI(x,y)])

static void Swap(double u[], int i, int o) {
  double t = u[i];
  u[i] = u[o];
  u[o] = t;
}

// x-component of normal, anti-gradient of u
static double NX(const double u[]) {
  const int w = SW;
  return GET(w - 1, w) - GET(w + 1, w);
}

// y-component of normal, anti-gradient of u
static double NY(const double u[]) {
  const int w = SW;
  return GET(w, w - 1) - GET(w, w + 1);
}

static void FlipX(double u[]) {
  for (int y = 0; y < SA; ++y) {
    for (int x = 0; x < SW; ++x) {
      Swap(u, GETI(x, y), GETI(SA - x - 1, y));
    }
  }
}

static void FlipY(double u[]) {
  for (int y = 0; y < SW; ++y) {
    for (int x = 0; x < SA; ++x) {
      Swap(u, GETI(x, y), GETI(x, SA - y - 1));
    }
  }
}

static void Trans(double u[]) {
  for (int y = 0; y < SA; ++y) {
    for (int x = 0; x < y; ++x) {
      Swap(u, GETI(x, y), GETI(y, x));
    }
  }
}

static double Mean(const double u[]) {
  double s = 0.;
  for (int i = 0; i < DIM; ++i) {
    s += u[i];
  }
  return s / DIM;
}

static void FlipVal(double u[]) {
  for (int i = 0; i < DIM; ++i) {
    u[i] = 1. - u[i];
  }
}

static double ApplySymm(double u[]) {
  int w = SW;

  double q = 1.;
  if (Mean(u) > 0.5) {
    FlipVal(u);
    q = -1.;
  }

  if (NX(u) < 0.) {
    FlipX(u);
  }

  if (NY(u) < 0.) {
    FlipY(u);
  }

  if (NY(u) < NX(u)) {
    Trans(u);
  }

#ifndef NDEBUG
  if (!(NX(u) >= 0. && NY(u) >= 0. && NX(u) <= NY(u))) {
    fprintf(stderr, "%s:%d nx=%g ny=%g\n", __FILE__, __LINE__, NX(u), NY(u));
    abort();
  }
#endif

  return q;
}

static double Curv(Point point, scalar c) {
  static double u[DIM];
  int i = 0;
  foreach_neighbor(SW)
  {
    u[i++] = c[];
  }
  double q = ApplySymm(u);
  return q * Eval(u) / Delta;
}


trace
cstats curvature_learn(struct Curvature p)
{
  scalar c = p.c, kappa = p.kappa;
  double sigma = p.sigma ? p.sigma : 1.;

#if TREE
  kappa.refine = kappa.prolongation = curvature_prolongation;
  kappa.restriction = curvature_restriction;
#endif

  scalar k[];
  scalar_clone (k, kappa);

  foreach() {
    if (!interfacial (point, c)) {
      k[] = nodata;
    } else  {
      k[] = Curv(point, c);
    }
  }

  boundary ({k});

  foreach () {
    double kf = k[];
    if (kf == nodata)
      kappa[] = nodata;
    else if (p.add)
      kappa[] += sigma*kf;
    else
      kappa[] = sigma*kf;
  }

  boundary ({kappa});

  return (cstats){0, 0, 0, 0};
}

trace
cstats curvature_orig(struct Curvature p) {
  return curvature(p);
}

#define curvature curvature_learn

