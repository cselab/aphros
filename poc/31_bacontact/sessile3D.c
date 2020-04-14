#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "vof.h"

scalar f[], * interfaces = {f};

#include "tension.h"

double theta0 = 45;

vector h[];
h.t[back] = contact_angle (theta0*pi/180.);
h.r[back] = contact_angle (theta0*pi/180.);

#define MAXLEVEL 5

int main()
{

  const face vector muc[] = {.1,.1,.1};
  mu = muc;
  f.height = h;
  f.sigma = 1.;
  N = 1 << MAXLEVEL;
  for (theta0 = 30; theta0 <= 150; theta0 += 30)
    run();
}

event init (t = 0)
{
  fraction (f, - (sq(x) + sq(y) + sq(z) - sq(0.5)));
}

event logfile (i += 10; t <= 10)
{
  scalar kappa[];
  cstats cs = curvature (f, kappa);
  foreach()
    if (f[] <= 1e-3 || f[] >= 1. - 1e-3)
      kappa[] = nodata;
  stats s = statsf (kappa);
  fprintf (fout, "%g %g %g %g %g %g %g %d %d %d %d\n", t, normf(u.x).max,
	   s.min, s.sum/s.volume, s.max, s.stddev, statsf(f).sum,
	   cs.h, cs.f, cs.a, cs.c);
  fflush (fout);
  if (s.stddev < 1e-2)
    return 1;
}

event end (t = end)
{
  scalar kappa[];
  curvature (f, kappa);
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 4.*statsf(f).sum;
  fprintf (ferr, "%d %g %.5g %.3g %.5g %.3g %.5g\n",
	   N, theta0, R, s.stddev, V, normf(u.x).max, t);
}

#if TREE
event adapt (i++) {
  scalar f1[];
  foreach()
    f1[] = f[];
  boundary ({f1});
  adapt_wavelet ({f1}, (double[]){1e-3}, minlevel = 3, maxlevel = MAXLEVEL);
}
#endif
