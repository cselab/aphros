#include "par.h"

#include ".u/curv/select.h"
#include ".u/io/io.h"
#include ".u/io/iompi.h"
#include ".u/bashape.h"
#include ".u/bah5.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "vof.h"
#include "tension.h"

int nxexp = REFINE;
int argnx;

#define ONROOT int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); if (rank == 0)

double sqr(double a) {
  return a * a;
}


WALLX
WALLY
#if dimension == 3
WALLZ
#endif

int main() {
  argnx = (1 << nxexp);

  init_grid(argnx);
  origin (0.,0.,0.);
  size(EXTENT);

  PERX
  PERY
#if dimension == 3
  PERZ
#endif
  
  rho1 = RHO2; 
  rho2 = RHO1; 

  mu1 = MU2;
  mu2 = MU1, 


  //f.height = h;
  f.sigma = SIGMA;

#ifdef CURV_PARTSTR
#ifdef PS_Np
  kPartstr.Np = PS_Np;
#endif
#ifdef PS_Ns
  #if dimension == 3
  kPartstr.Ns = PS_Ns;
  #endif
#endif
#ifdef PS_Hp
  kPartstr.Hp = PS_Hp;
#endif
#ifdef PS_eps
  kPartstr.eps = PS_eps;
#endif
#ifdef PS_itermax
  kPartstr.itermax = PS_itermax;
#endif
#ifdef PS_eta
  kPartstr.eta = PS_eta;
#endif
#ifdef PS_circ
  kPartstr.circ = PS_circ;
#endif
#endif

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

  const int MAX_SHAPES = 10000;
  int nb = MAX_SHAPES;
  Shape bb[MAX_SHAPES];
  nb = ReadList("../b.dat", bb, nb);

  scalar ft[];
  for (int i = 0; i < nb; ++i) {
    Shape b = bb[i];
#if dimension == 2
    b.z = 0;
    b.u = 0;
    b.v = 0;
    b.w = 1;
    b.type = SHAPE_C;
#endif
    CreateField(b, ft);
    foreach () {
      f[] = max(f[], ft[]);
    }
  }
  boundary({f});
}

event out (t += DUMPDT ; t <= TMAX + DUMPDT) {
  ONROOT fprintf(stderr, "dump i=%05d t=%g dt=%g \n", i, t ,dt);

  static int frame = 0;
  //scalar * a = {u, p, f};
  scalar * a = {f};

  char name[1000];

#if dimension == 2
  sprintf(name, "o/%d/u_%04d.vtk", pid(), frame);
  FILE * fp = fopen(name, "w");
  io(a, fp);
  fclose(fp);
#elif dimension == 3
  sprintf(name, "%%s_%04d", frame);
  scalar vf[];
  foreach () {
    vf[] = f[];
  }
  bah5_list({vf}, name); /* '%s is replaced by a field name */
#endif

  ++frame;
}

event statout (i += 10 ; t <= TMAX) {
  ONROOT fprintf(stderr, "step i=%05d t=%g dt=%g \n", i, t ,dt);
}

event logfile (i += 1 ; t <= TMAX) {
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

  sb = fmax(sb, 1e-10);

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
