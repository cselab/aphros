#include "fractions.h"
#include "curvature.h"

#include <unistd.h>

#ifndef POINTXYZ
#define POINTXYZ
#endif

#ifndef foreach_end
#define foreach_end
#endif

#ifndef foreach_neighbor_end
#define foreach_neighbor_end
#endif

#ifndef CELLIDX
#define CELLIDX
#endif

#define SA (SW*2+1)


#include "learn.h"

#define GET(i,j) (u[(i)*SA+(j)])

static void Swap(double u[], int i, int o) {
  double t = u[i];
  u[i] = u[o];
  u[o] = t;
}

static double A(const double u[]) {
  const int w = SW;
  double nx = GET(w, w + 1) - GET(w, w - 1);
  double ny = GET(w + 1, w) - GET(w - 1, w);
  return atan2(ny, nx) * 180. / M_PI;
}

static void FlipX(double u[]) {
  for (int i = 0; i < SW; ++i) {
    for (int j = 0; j < SA; ++j) {
      Swap(u, j * SA + i, j * SA + (SA - i - 1));
    }
  }
}

static void FlipY(double u[]) {
  for (int i = 0; i < SA; ++i) {
    for (int j = 0; j < SW; ++j) {
      Swap(u, j * SA + i, (SA - j - 1) * SA + i);
    }
  }
}

static void Trans(double u[]) {
  for (int i = 0; i < SA; ++i) {
    for (int j = 0; j < i; ++j) {
      Swap(u, j * SA + i, i * SA + j);
    }
  }
}

static void ApplySymm(double u[]) {
  int w = SW;

  if (A(u) < 0) {
    FlipY(u);
  }

  if (A(u) > 90) {
    FlipX(u);
  }

  if (A(u) > 45) {
    Trans(u);
  }

  double a = A(u);
  const double TH = 1e-14;
  if (!(a >= -TH && a <= 45 + TH)) {
    fprintf(stderr, "%s:%d %g\n", __FILE__, __LINE__, a);
    abort();
  }
}

static double partstr(Point point, scalar c) {
  static double u[DIM];
  int i = 0;
  foreach_neighbor(SW)
  {
    u[i++] = c[];
  }
  ApplySymm(u);
  return Eval(u) / Delta;
}



#ifndef NOBA
trace
cstats curvature_learn(struct Curvature p)
{
  /*
  {
    scalar kappa = p.kappa;
    foreach () {
      kappa[CELLIDX] = 0;
    }
    return (cstats){0,0,0,0};
  }
  */
  scalar c = p.c, kappa = p.kappa;
  double sigma = p.sigma ? p.sigma : 1.;
  int sh = 0, sf = 0, sa = 0, sc = 0;
  vector ch = c.height, h = automatic (ch);
  if (!ch.x.i)
    heights (c, h);

#if TREE
  kappa.refine = kappa.prolongation = curvature_prolongation;
  kappa.restriction = curvature_restriction;
#endif

  scalar k[];
  scalar_clone (k, kappa);

  foreach(reduction(+:sh) reduction(+:sc)) {
    if (!interfacial (point, c)) {
      k[CELLIDX] = nodata;
    } else  {
      k[CELLIDX] = partstr(point, c);
      sc++;
    }
  }
  foreach_end
  boundary ({k});

  foreach () {
    double kf = k[CELLIDX];
    if (kf == nodata)
      kappa[CELLIDX] = nodata;
    else if (p.add)
      kappa[CELLIDX] += sigma*kf;
    else
      kappa[CELLIDX] = sigma*kf;
  }
  foreach_end
  boundary ({kappa});

  return (cstats){sh, sf, sa, sc};
}

trace
cstats curvature_orig(struct Curvature p) {
  return curvature(p);
}

#define curvature curvature_learn

#endif

