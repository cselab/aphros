#include "grid/multigrid3D.h"
#include "fractions.h"
#include "curvature_select.h"
#include "io/io.h"
#include <vof.h>

#include <assert.h>
#define myassert(EX) (void)((EX) || (__assert (#EX, __FILE__, __LINE__),0))

#include <mpi.h>

#define ONROOT int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); if (rank == 0)

typedef struct Shape Shape;
struct Shape {
  double x, y, z;
  double rx, ry, rz;
  double u, v, w;
  int type;
};
struct Shape b;

enum {BUB_S, BUB_C, BUB_E};
int nxexp;
int argnx;

void CreateField(Shape b, scalar c) {
  foreach() {
    double h = Delta;
    c[] = vof_cylinder((b.x - x)/h , (b.y - y)/h, (b.z - z)/h,
                       b.rx/h, b.u, b.v, b.w);
  }
}

int Inside(Shape b, double x, double y, double z, double h) {
  return sq(x -  b.x) + sq(y -  b.y) + sq(z -  b.z) < sq(b.rx + h);
}

int Read(Shape b, FILE *f) {
  fscanf(f, "%lf %lf %lf %lf", &b.x, &b.y, &b.z, &b.rx);
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
  Read(b, fb);
  fclose(fb);
  origin (0.,0.,0.);
  scalar vf[];
  CreateField(b, vf);
  scalar k[];
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
        if (Inside(b, x, y, z, Delta))
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
