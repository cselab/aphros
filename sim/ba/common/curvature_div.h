#include "fractions.h"
#include "curvature.h"

#include <unistd.h>

#define SA 3
#define SW 1
#define DIM (SW*SW)
#define GET(x,y) (u[(y)*SA+(x)])
#define U(x,y) (u[(y+SW)*SA+(x+SW)])

static double Norm(double x, double y) {
  return sqrt(x * x + y * y);
}

static double Curv(Point point, scalar c) {
  static double u[DIM];
  int i = 0;
  foreach_neighbor(SW)
  {
    u[i++] = c[];
  }

  // gradient
  double gxpx = U(1,0) - U(0,0);
  double gxpy = ((U(0,1) + U(1,1)) - (U(0,-1) + U(1,-1))) * 0.25;

  double gxmx = U(0,0) - U(-1,0);
  double gxmy = ((U(-1,1) + U(0,1)) - (U(-1,-1) + U(0,-1))) * 0.25;

  double gypx = ((U(1,0) + U(1,1)) - (U(-1,0) + U(-1,1))) * 0.25;
  double gypy = U(0,1) - U(0,0);

  double gymx = ((U(1,-1) + U(1,0)) - (U(-1,-1) + U(-1,0))) * 0.25;
  double gymy = U(0,0) - U(0,-1);

  // length
  double gxp = Norm(gxpx, gxpy);
  double gxm = Norm(gxmx, gxmy);
  double gyp = Norm(gypx, gypy);
  double gym = Norm(gymx, gymy);

  // divergence of unit normal
  double d = gxpx / gxp - gxmx / gxm + gypy / gyp - gymy / gym;

  return d / Delta;
}


trace
cstats curvature_div(struct Curvature p)
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
    if (!interfacial(point, c)) {
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

#define curvature curvature_div

