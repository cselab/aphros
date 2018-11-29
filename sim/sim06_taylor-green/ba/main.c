#include "par.h"

#define FILTERED

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "vtk.h"
#include <mpi.h>

#include "io/iompi.h"

scalar f0[];
double uin = INV;
u.n[top]  = dirichlet(f0[] * uin);
f[top]    = 1. - f0[];
//u.t[top]  = dirichlet(0);
//u.r[top]  = dirichlet(0);
//p[top]    = neumann(0);


u.n[left]  = (y > OUTY0 && y < OUTY1 ? neumann(0) : dirichlet(0));
u.n[right]  = (y > OUTY0 && y < OUTY1 ? neumann(0) : dirichlet(0));
u.n[front]  = (y > OUTY0 && y < OUTY1 ? neumann(0) : dirichlet(0));
u.n[back]  = (y > OUTY0 && y < OUTY1 ? neumann(0) : dirichlet(0));



/*
u.n[left]  = neumann(0);
u.n[right]  = neumann(0);
u.n[front]  = neumann(0);
u.n[back]  = neumann(0);
*/

//u.n[bottom] = neumann(0);
p[bottom]   = dirichlet(0);


int main() {
  init_grid(NX);

  origin(0., 0., 0.);
  dimensions(nx = 4, ny = 16, nz = 4);
  //dimensions(nx = 2, ny = 8, nz = 2);
  //size(1.);

  rho1 = RHO1; 
  rho2 = RHO2; 

  mu1 = MU1;
  mu2 = MU2, 

  f.sigma = SIGMA;

  //f.sigma = 0;
  //rho2 = RHO1; 
  //mu2 = MU1, 


  run();
}

event init (i = 0) {
  double x0 = INX0;
  double x1 = INX1;
  double k = sqrt(4. / pi);
  double R = (x1 - x0) * 0.5 * k;
  //fraction(f0, -min(min(min(x - x0, x1 - x), z - x0), x1 - z));
  fraction(f0, sq(R) - (sq(x - 0.5) + sq(z - 0.5)));
  fraction(f, y - (EXTENT - AIRGAP));

  DT = 1e-4
}


/*
// gravity
event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] += GY;
}
*/

event updatedt (i++) {
  if (t < 1e-2) {
    DT = 1e-4;
    //TOLERANCE = 1e-2;
  } else {
    DT = 1e-3 * 0.2;
    //TOLERANCE = 1e-2;
  }
}


event out (t += DUMPDT ; t <= TMAX) {
  static int frame = 0;
  char name[1000];
  sprintf(name, "o/%d/u_%04d.vtk", pid(), frame);
  ++frame;
  //scalar * a = {u, p, f};
  scalar * a = {f};
  iompi(a, name);
}

event statout (i += 1) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
  printf("step i=%05d t=%g dt=%g \n", i, t ,dt);
  }
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
    vlmx = fmax(vlmx, fabs(u.x[]) * f[]);
    vlmy = fmax(vlmy, fabs(u.y[]) * f[]);
    vlmz = fmax(vlmz, fabs(u.z[]) * f[]);

    double dv = f[]*dv();
    vl2x += sq(u.x[]) * dv;
    vl2y += sq(u.y[]) * dv;
    vl2z += sq(u.z[]) * dv;
  }
  static char* fn = "o/sc";
  static FILE* f = NULL;

  sb = fmax(sb, 1e-10);

  //int rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //if (rank == 0) {

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
  //}
}
