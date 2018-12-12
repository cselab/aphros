#include "par.h"

#include "fractions.h"
#include "curvature.h"

#include "io/io.h"

double sqr(double a) {
  return a * a;
}

struct Bub {
  double x;
  double y;
  double z;
  double r;
} b;

double ifr2(double x, double y) {
  double r2 = sq(x - b.x) + sq(y - b.y);
  return sq(b.r) - r2;
}

double ifr3(double x, double y, double z) {
  double r2 = sq(x - b.x) + sq(y - b.y) + sq(z - b.z);
  return sq(b.r) - r2;
}

double ls(double x, double y, double z) {
  return ifr3(x, y, z);
}

void fraction2(scalar c) {
  foreach_cell() {
    double h = Delta;
    double hh = h * 0.5;
    if (
        ls(x, y, z) * ls(x + h, y + h, z + h) > 0. && 
        ls(x, y, z) * ls(x + h, y + h, z - h) > 0. && 
        ls(x, y, z) * ls(x + h, y - h, z + h) > 0. && 
        ls(x, y, z) * ls(x + h, y - h, z - h) > 0. && 
        ls(x, y, z) * ls(x - h, y + h, z + h) > 0. && 
        ls(x, y, z) * ls(x - h, y + h, z - h) > 0. && 
        ls(x, y, z) * ls(x - h, y - h, z + h) > 0. && 
        ls(x, y, z) * ls(x - h, y - h, z - h) > 0.
        ) {
      c[] = (ls(x, y, z) > 0. ? 1. : 0.);
    } else {
      double dx = 1e-3 * h;
      double q = dx * 0.5;
      double nx = (ls(x + q, y, z) - ls(x - q, y, z)) / dx;
      double ny = (ls(x, y + q, z) - ls(x, y - q, z)) / dx;
      double nz = (ls(x, y, z + q) - ls(x, y, z - q)) / dx;
      double a = ls(x, y, z);
      c[] = rectangle_fraction(
        (coord){nx, ny, nz}, a, 
        (coord){-hh, -hh, -hh}, 
        (coord){hh, hh, hh}); 
    }
  }

  double a = 0.;
  foreach() {
    double h = Delta;
    a += c[] * h * h * h;
  }
  printf("%.16g", a);
}


int main() {
  init_grid(1 << REFINE);

  FILE* fb = fopen("../b.dat", "r");
  fscanf(fb, "%lf %lf %lf %lf", &b.x, &b.y, &b.z, &b.r);
  fclose(fb);

  origin (0.,0.,0.);

  scalar vf[]; // volume fraction

  /*
#if dimension == 2
  fraction(vf, ifr2(x, y));
#elif dimension == 3
  fraction(vf, ifr3(x, y, z));
#endif
*/

  fraction2(vf);

  scalar k[]; // curvature
  //curvature(vf, k);
  
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

