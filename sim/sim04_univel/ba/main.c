#include "par.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "contact.h"
#include "vof.h"
#include "tension.h"

#include "io/iompi.h"
#include "io/io.h"

double sqr(double a) {
  return a * a;
}

double ifr2(double x, double y) {
  double r = sq(x - BCX) + sq(y - BCY);
  double r2 = sq(x - BC2X) + sq(y - BC2Y);
  return fmax(sq(BR) - r, sq(BR2) - r2);
}

double ifr3(double x, double y, double z) {
  double r = sq(x - BCX) + sq(y - BCY) + sq(z - BCZ);
  double r2 = sq(x - BC2X) + sq(y - BC2Y) + sq(z - BC2Z);
  return fmax(sq(BR) - r, sq(BR2) - r2);
}

vector h[];

WALLX
WALLY
WALLZ

int main() {
  init_grid(1 << REFINE);

  origin (0.,0.,0.);

  PERX
  PERY
  PERZ
  
  rho1 = RHO2; 
  rho2 = RHO1; 

  mu1 = MU2;
  mu2 = MU1, 

  f.height = h;
  f.sigma = SIGMA;

  run();
}

event init (i = 0) {
  foreach() {
    u.x[] = VELX;
    u.y[] = VELY;
#if dimension == 3
    u.z[] = VELZ;
#endif
  }

#if dimension == 2
  fraction(f, ifr2(x, y));
#elif dimension == 3
  fraction(f, ifr3(x, y, z));
#endif
}

event out (t += DUMPDT ; t <= TMAX) {
  static int frame = 0;
  ++frame;
  //scalar * a = {u, p, f};
  scalar * a = {f};

  char name[1000];
  sprintf(name, "o/%d/u_%04d.vtk", pid(), frame);

#if dimension == 2
  FILE * fp = fopen(name, "w");
  io(a, fp);
  fclose(fp);
#elif dimension == 3
  iompi(a, name);
#endif
}

event statout (i += 10) {
  printf("step i=%d\n", i);
}

event logfile (i += 1) {
  double xb = 0., yb = 0., zb = 0., sb = 0.;
  double vbx = 0., vby = 0., vbz = 0.;
  double p0 = 1e10, p1 = -1e10;
  foreach(reduction(+:xb) reduction(+:yb) reduction(+:zb)
   reduction(+:vbx) reduction(+:vby) reduction(+:vbz)
   reduction(+:sb) 
   reduction(min:p0) reduction(max:p1)
   ) {
    double dv = f[]*dv();
    xb += x*dv;
    yb += y*dv;
    zb += z*dv;
    vbx += u.x[]*dv;
    vby += u.y[]*dv;
    vbz += u.z[]*dv;
    sb += dv;
    p0 = fmin(p0, p[]);
    p1 = fmax(p1, p[]);
  }

  double vlmx = 0., vlmy = 0., vlmz = 0.;
  double vl2x = 0., vl2y = 0., vl2z = 0.;
  foreach(reduction(+:vl2x) reduction(+:vl2y) reduction(+:vl2z)
   reduction(max:vlmx) reduction(max:vlmy) reduction(max:vlmz)
   ) {
    vlmx = fmax(vlmx, fabs(u.x[] - VELX) * f[]);
    vlmy = fmax(vlmy, fabs(u.y[] - VELY) * f[]);
    vlmz = fmax(vlmz, fabs(u.z[] - VELZ) * f[]);

    double dv = f[]*dv();
    vl2x += sqr(u.x[] - VELX) * dv;
    vl2y += sqr(u.y[] - VELY) * dv;
    vl2z += sqr(u.z[] - VELZ) * dv;
  }
  static char* fn = "o/sc";
  static FILE* f = NULL;

  if (!f) {
    f = fopen(fn, "w");
    fprintf(f, "t m2 c2x c2y c2z v2x v2y v2z p0 p1 pd vlmx vlmy vlmz vl2x vl2y vl2z\n");
  } 

  fprintf(f,
    "%.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f\n",
    t, sb,
    xb/sb, yb/sb, zb/sb,
    vbx/sb, vby/sb, vbz/sb, p0, p1, p1 - p0,
    vlmx, vlmy, vlmz, 
    sqrt(vl2x/sb), sqrt(vl2y/sb), sqrt(vl2z/sb)
    );

  fflush(f);
}
