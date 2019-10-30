#include <mpi.h>

#include "par.h"

#include ".u/bashape.h"
#include ".u/bah5.h"
#include ".u/curv/select.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "vof.h"
#include "tension.h"

#ifndef BRZ
#define BRZ BR
#endif

#ifndef BRY
#define BRY BR
#endif

#define EXTENTY (EXTENT * 0.5)

#define ONROOT int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); if (rank == 0)

double ifr3(double x, double y, double z) {
  double r = sq((x - BCX) / BR) + sq((y - BCY) / BRY) + sq((z - BCZ) / BRZ);
  return sq(BR) * (1 - r);
}

double Vel(double y) {
  double dy = (y - BCY) / (2 * BR);
  return VELC + SHEAR * dy;
}

u.r[bottom] = dirichlet(Vel(y));
u.t[bottom] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.r[top] = dirichlet(Vel(y));
u.t[top]  = dirichlet(0.);
u.n[top]  = dirichlet(0.);

#define P(X) { ONROOT fprintf(stderr, #X "=%g\n", (double)X); }

int main() {
  init_grid(NX);

  P(RHO1)
  P(RHO2)
  P(MU1)
  P(MU2)
  P(SIGMA)
  P(DUMPDT)
  P(TMAX)
  P(EXTENT)
  P(EXTENTY)
  P(WE)
  P(RE)
  P(BCX)
  P(BCY)
  P(BCZ)
  P(BR)
  P(SHEAR)
  P(VEL)
  P(VELC)

  periodic(right);
  periodic(front);

  origin(0., 0., 0.);
  MPIDIM

  size(EXTENT);

  rho1 = RHO2;
  rho2 = RHO1;

  mu1 = MU2;
  mu2 = MU1;

  f.sigma = SIGMA;

  run();
}

event init (i = 0) {
  fraction(f, ifr3(x, y, z));
  foreach() {
    //u.x[] = (Vel(y)) * (1. - f[]);
  }
}


event out (t += DUMPDT ; t <= TMAX) {
#if !NOIO
  static int frame = 0;
  char name[1000];
  scalar * a = {u, f};
  //scalar * a = {f};
  sprintf(name, "%%s_%04d", frame);
  bah5_list(a, name);
  ++frame;
#endif
  ONROOT fprintf(stderr, "dump %s step i=%05d t=%g\n", name, i, t);
}

event statout (i += 1) {
  ONROOT fprintf(stderr, "step i=%05d t=%g dt=%g \n", i, t ,dt);
}

event logfile (i += 1 ; t <= TMAX) {
  // center of mass and volume
  double xb = 0, yb = 0, zb = 0, sb = 0;
  // gyration tensor
  double xxb = 0, xyb = 0, xzb = 0, yyb = 0, yzb = 0, zzb = 0;
  // mean velocity
  double vbx = 0, vby = 0, vbz = 0;
  // min and max pressure
  double p0 = 1e10, p1 = -1e10;
  foreach(
   reduction(+:xb) reduction(+:yb) reduction(+:zb)
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
  foreach(
   reduction(+:xxb) reduction(+:xyb) reduction(+:xzb)
   reduction(+:yyb) reduction(+:yzb) reduction(+:zzb)
   ) {
    double dv = f[]*dv();
    xxb += x*x*dv;
    xyb += x*y*dv;
    xzb += x*z*dv;
    yyb += y*y*dv;
    yzb += y*z*dv;
    zzb += z*z*dv;
  }

  sb = fmax(sb, 1e-10);

  static char* fn = "o/sc";
  static FILE* f = NULL;

  if (!f) {
    f = fopen(fn, "w");
    fprintf(f, "t m2 c2x c2y c2z v2x v2y v2z p0 p1 pd xx xy xz yy yz zz\n");
  } 

  fprintf(f,
    "%.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f\n",
    t, sb,
    xb/sb, yb/sb, zb/sb,
    vbx/sb, vby/sb, vbz/sb, p0, p1, p1 - p0,
    xxb/sb, xyb/sb, xzb/sb, yyb/sb, yzb/sb, zzb/sb
    );

  fflush(f);
}
