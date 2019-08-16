#include "fractions.h"
#include "curvature.h"

#include <unistd.h>

#include "learn_eval.h"


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
  const int sdim = SX * SY * SZ;
  double s = 0.;
  for (int i = 0; i < sdim; ++i) {
    s += u[i];
  }
  return s / sdim;
}

static void FlipVal(double u[]) {
  const int sdim = SX * SY * SZ;
  for (int i = 0; i < sdim; ++i) {
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

static double ApplySymm2(double u[]) {
#if dimension == 2
  return ApplySymm(u);
#else
  static const int sw = SW;
  static const int sa = SW * 2 + 1;
  double v[sa * sa * sa];

  {
    int j = 0;
    for (int z = 0; z < sa; ++z) {
      int i = 0;
      for (int y = 0; y < SY; ++y) {
        for (int x = 0; x < SX; ++x) {
          v[j++] = u[i++];
        }
      }
    }
  }

  double q = ApplySymm(v);
  {
    int i = 0;
    int j = 0;
    for (int y = 0; y < SY; ++y) {
      for (int x = 0; x < SX; ++x) {
        u[i++] = v[j++];
      }
    }
  }
  return q;
#endif
}

// u: volume fraction in 3D stencil
// uu: volume fraction in 2D stencil
// d0,d1: slice directions (0,1,2)
static void Slice(const double u[], double uu[], int d0, int d1) {
  static const int sw = SW;
  static const int sa = SW * 2 + 1;
  int ii = 0;
  for (int i1 = 0; i1 < sa; ++i1) {
    for (int i0 = 0; i0 < sa; ++i0) {
      int q[3] = {sw, sw, sw}; // central cell
      q[d0] = i0;
      q[d1] = i1;
      int i = sa * sa * q[2] + sa * q[1] + q[0];
      uu[ii] = u[i];
      ++ii;
    }
  }
  assert(ii == sa * sa);
}


static double EvalSlice(const double u[], int d0, int d1) {
  double uu[SX * SX];
  Slice(u, uu, d0, d1);
  double q = ApplySymm2(uu);
  return q * Eval(uu);
}


static double Curv(Point point, scalar c) {
  double u[SX * SX * SX];
  int i = 0;
  foreach_neighbor(SW) {
    u[i++] = c[];
  }
  double q = ApplySymm(u);

#if DIM == SX * SY * SZ
  double k = Eval(u);
#elif DIM == SX * SY
  double k = EvalSlice(u, 0, 1) + EvalSlice(u, 0, 2);
#else
  #error Incorrect input size DIM and SX*SY*SZ
#endif
  return q * k / Delta;
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

