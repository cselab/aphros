#include "fractions.h"
#include "curvature.h"

#include "io/io.h"

#include <assert.h>
#define myassert(EX) (void)((EX) || (__assert (#EX, __FILE__, __LINE__),0))

#include <mpi.h>

#define ONROOT int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); if (rank == 0)

double sqr(double a) {
  return a * a;
}

struct Bub {
  double x;
  double y;
  double z;
  double r;
} b;

int nxexp;
int argnx;

double ifr2(double x, double y) {
  double r2 = sq(x - b.x) + sq(y - b.y);
  return sq(b.r) - r2;
}

double ifr3(double x, double y, double z) {
  double r2 = sq(x - b.x) + sq(y - b.y) + sq(z - b.z);
  return sq(b.r) - r2;
}

double ls(double x, double y, double z) {
#if dimension == 2
  return ifr2(x, y);
#elif dimension == 3
  return ifr3(x, y, z);
#endif
}

void fraction2(scalar c) {
  foreach() {
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

}

void ReadField(scalar c, char* fn) {
  FILE* f = fopen(fn, "r");

  int nx, ny, nz;
  fscanf(f, "%d %d %d", &nx, &ny, &nz);

  myassert(nx == argnx);
  myassert(ny == argnx);

#if dimension == 2
  myassert(nz == 1);
#elif dimension == 3
  myassert(nz == argnx);
#endif

  double* uu = malloc(sizeof(double) * nx * ny * nz);

  for (int z = 0; z < nz; ++z) {
    for (int y = 0; y < ny; ++y) {
      for (int x = 0; x < nx; ++x) {
        double a;
        fscanf(f, "%lf", &a);
        uu[z * ny * nx + y * nx + x] = a;
      }
    }
  }
  fclose(f);

  int i = 0;
  foreach() {
    double h = Delta;
    double hmin = 1. / argnx;
    int ix = max(0, min(x / hmin, nx - 1));
    int iy = max(0, min(y / hmin, ny - 1));
    int iz = max(0, min(z / hmin, nz - 1));
    c[] = uu[iz * ny * nx + iy * nx + ix];
  }

  free(uu);

  boundary ({c});
}


int main() {
  {
    FILE* q = fopen("nxexp", "r");
    fscanf(q, "%d", &nxexp);
    fclose(q);
  }

  argnx = (1 << nxexp);

  init_grid(argnx);

  FILE* fb = fopen("b.dat", "r");
  fscanf(fb, "%lf %lf %lf %lf", &b.x, &b.y, &b.z, &b.r);
  fclose(fb);

  origin (0.,0.,0.);

  scalar vf[]; // volume fraction

  ReadField(vf, "../ch/vf_0000.dat");

  scalar k[]; // curvature
  curvature(vf, k);

#if dimension == 2
  double kc = 1.;
#elif dimension == 3
  double kc = 0.5;
#endif
  
  {
    FILE* q = fopen("ok", "w");
    foreach() {
      if (vf[] > 0. && vf[] < 1.) {
        fprintf(q, "%.16g\n", k[] * kc);
      }
    }
    fclose(q);
  }

  {
    FILE* q = fopen("ovf", "w");
    double s = 0.;
    foreach() {
      double h = Delta;
      s += vf[] * h * h * h;
    }
    fprintf(q, "%.16g\n", s);
    fclose(q);
  }


  {
    FILE* q = fopen("u.vtk", "w");
    io({vf, k}, q);
    fclose(q);
  }
}

