#include "par.h"

#include "fractions.h"
#include "curvature.h"

#include "io/io.h"

double sqr(double a) {
  return a * a;
}

double ifr2(double x, double y, 
    double bx, double by, double br) {
  double r = sq(x - bx) + sq(y - by);
  return sq(br) - r;
}

double ifr3(double x, double y, double z, 
    double bx, double by, double bz, double br) {
  double r = sq(x - bx) + sq(y - by) + sq(z - bz);
  return sq(br) - r;
}


int main() {
  init_grid(1 << REFINE);

  double bx, by, bz, br;
  FILE* fb = fopen("../b.dat", "r");
  fscanf(fb, "%lf %lf %lf %lf", &bx, &by, &bz, &br);
  fclose(fb);

  origin (0.,0.,0.);

  scalar vf[]; // volume fraction

#if dimension == 2
  fraction(vf, ifr2(x, y, bx, by, br));
#elif dimension == 3
  fraction(vf, ifr3(x, y, z, bx, by, bz ,br));
#endif

  scalar k[]; // curvature
  curvature(vf, k);
  
  {
    FILE* q = fopen("t", "w");
    foreach() {
      if (vf[] > 0. && vf[] < 1.) {
        fprintf(q, "%.16g\n", k[]*0.5);
      }
    }
    fclose(q);
  }

  {
    FILE* q = fopen("u.vtk", "w");
    io({vf, k}, q);
    fclose(q);
  }
}

