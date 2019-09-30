#include ".u/curv/select.h"
#include ".u/io/io.h"
#include ".u/bashape.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "vof.h"

#define DUMPDT 0.01
#define TMAX 1

int main() {
  init_grid(64);
  origin (0.,0.,0.);

  rho1 = 1;
  rho2 = 0.01;

  mu1 = 0;
  mu2 = 0;

  run();
}

event init (i = 0) {
  foreach() {
    u.x[] = 1;
    u.y[] = 0;
  }
}

event out (t += DUMPDT ; t <= TMAX + DUMPDT) {
  fprintf(stderr, "dump i=%05d t=%g dt=%g \n", i, t ,dt);

  static int frame = 0;
  scalar * a = {u, p, f};

  char name[1000];

  sprintf(name, "o/%d/u_%04d.vtk", pid(), frame);
  FILE * fp = fopen(name, "w");
  io(a, fp);
  fclose(fp);

  ++frame;
}

event statout (i += 10 ; t <= TMAX) {
  fprintf(stderr, "step i=%05d t=%g dt=%g \n", i, t ,dt);
}
