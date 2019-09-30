#include ".u/curv/select.h"
#include ".u/io/io.h"
#include ".u/io/iompi.h"
#include ".u/bashape.h"
#include ".u/bah5.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "vof.h"
#include "tension.h"

int main() {
  init_grid(64);
  origin (0.,0.,0.);

  rho1 = RHO2; 
  rho2 = RHO1; 

  mu1 = MU2;
  mu2 = MU1, 

  f.sigma = SIGMA;
  run();
}

event init (i = 0) {
  foreach() {
    u.x[] = 1;
    u.y[] = 0;
  }
}

event out (t += DUMPDT ; t <= TMAX + DUMPDT) {
  ONROOT fprintf(stderr, "dump i=%05d t=%g dt=%g \n", i, t ,dt);

  static int frame = 0;
  scalar * a = {u, p, f};

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
