#include "par.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "vtk.h"

#include "io/io.h"

int main() {
  init_grid(1 << REFINE);

  printf("dim=%d", dimension);

  origin (0.,0.,0.);
  foreach_dimension() {
    periodic (right);
  }
  
  rho1 = RHO1, mu1 = MU1;
  rho2 = RHO2, mu2 = RHO2, 
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
  fraction (f, sq(BR) - (sq(x - BCX) + sq(y - BCY)));
#elif dimension == 3
  fraction (f, sq(BR) - (sq(x - BCX) + sq(y - BCY) + sq(z - BCZ)));
#endif
}

event out (t += DUMPDT ; t < TMAX) {
  static int frame = 0;
  char name[1000];
  sprintf(name, "o/%d/u_%04d.vtk", pid(), frame);
  ++frame;
  FILE * fp = fopen(name, "w");
  scalar * a = {u, p, f};
  io(a, fp);
}
