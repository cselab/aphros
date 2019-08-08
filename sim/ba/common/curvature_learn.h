#include "fractions.h"
#include "curvature.h"

#include <unistd.h>

#include "learn.h"


#define SX (SW*2+1)
#define SY (SX)

#if dimension == 2
#define SZ 1
#define GETI(x,y,z) ((y)*SX + (x))
#else 
#define SZ (SX)
#define GETI(x,y,z) ((z)*SX*SY + (y)*SX + (x))
#endif

#define GET(x,y,z) (u[GETI(x,y,z)])

static void Swap(double u[], int i, int o) {
  double t = u[i];
  u[i] = u[o];
  u[o] = t;
}

// components of normal, anti-gradient of u
static double NX(const double u[]) {
  const int x = SX / 2;
  const int y = SY / 2;
  const int z = SZ / 2;
  return GET(x - 1, y, z) - GET(x + 1, y, z);
}
static double NY(const double u[]) {
  const int x = SX / 2;
  const int y = SY / 2;
  const int z = SZ / 2;
  return GET(x, y - 1, z) - GET(x, y + 1, z);
}
// returns 0 in 2D
static double NZ(const double u[]) {
  const int x = SX / 2;
  const int y = SY / 2;
  const int z = SZ / 2;
  return GET(x, y, z - 1) - GET(x, y, z + 1);
}

static void FlipX(double u[]) {
  for (int z = 0; z < SZ; ++z) {
    for (int y = 0; y < SY; ++y) {
      for (int x = 0; x < SX / 2; ++x) {
        Swap(u, GETI(x, y, z), GETI(SX - x - 1, y, z));
      }
    }
  }
}

static void FlipY(double u[]) {
  for (int z = 0; z < SZ; ++z) {
    for (int y = 0; y < SY / 2; ++y) {
      for (int x = 0; x < SX; ++x) {
        Swap(u, GETI(x, y, z), GETI(x, SY - y - 1, z));
      }
    }
  }
}

static void FlipZ(double u[]) {
  for (int z = 0; z < SZ / 2; ++z) {
    for (int y = 0; y < SY; ++y) {
      for (int x = 0; x < SX; ++x) {
        Swap(u, GETI(x, y, z), GETI(x, y, SZ - z - 1));
      }
    }
  }
}

static void TransXY(double u[]) {
  for (int z = 0; z < SZ; ++z) {
    for (int y = 0; y < SY; ++y) {
      for (int x = 0; x < y; ++x) {
        Swap(u, GETI(x, y, z), GETI(y, x, z));
      }
    }
  }
}

static void TransYZ(double u[]) {
  for (int z = 0; z < SZ; ++z) {
    for (int y = 0; y < z; ++y) {
      for (int x = 0; x < SX; ++x) {
        Swap(u, GETI(x, y, z), GETI(x, z, y));
      }
    }
  }
}

static void TransZX(double u[]) {
  for (int z = 0; z < SZ; ++z) {
    for (int y = 0; y < SY; ++y) {
      for (int x = 0; x < z; ++x) {
        Swap(u, GETI(x, y, z), GETI(z, y, x));
      }
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

  // make components positive
  if (NX(u) < 0.) FlipX(u);
  if (NY(u) < 0.) FlipY(u);
  if (NZ(u) < 0.) FlipZ(u);

  // order components as 0 <= nz <= ny <= nx
  if (NX(u) < NY(u)) TransXY(u);
  if (NY(u) < NZ(u)) TransYZ(u);
  if (NX(u) < NY(u)) TransXY(u);

#ifndef NDEBUG
  if (!(0. <= NZ(u) && NZ(u) <= NY(u) && NY(u) <= NX(u))) {
    fprintf(stderr, "%s:%d nx=%g ny=%g nz=%g\n", 
        __FILE__, __LINE__, NX(u), NY(u), NZ(u));
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

