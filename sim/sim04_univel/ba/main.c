#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "vtk.h"

#include "io/io.h"

#include "par.h"

int main() {
  init_grid(1 << REFINE);

  origin (0.,0.,0.);
  foreach_dimension() {
    periodic (right);
  }
  
  rho1 = RHO1, mu1 = MU1;
  rho2 = RHO2, mu2 = RHO2, 
  f.sigma = SIGMA;

  /**
  We reduce the tolerance on the Poisson and viscous solvers to
  improve the accuracy. */
  
  //DT = .1;

  run();
}

event init (i = 0) {
  foreach() {
    u.x[] = VELX;
    u.y[] = VELY;
  }

/*
  boundary ({p});
  foreach() {
    foreach_dimension() {
      g.x[] = - (p[1] - p[-1])/(2.*Delta);
    }
  }
  boundary ((scalar *){g});
*/

  fraction (f, sq(BR) - (sq(x - BCX) + sq(y - BCY) + sq(z - BCZ)));
}

event out (t += DUMPDT ; t < TMAX) {
  static int frame = 0;
  char name[1000];
  sprintf(name, "o/%d/u_%04d.vtk", pid(), frame);
  //sprintf(name, "%04d.vtk", frame);
  ++frame;
  FILE * fp = fopen(name, "w");
  scalar * a = {u, p, f};
  io(a, fp);
}
