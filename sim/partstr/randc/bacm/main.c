#include <assert.h>
#include "grid/multigrid3D.h"
#include "fractions.h"
#include "curvature_select.h"
#include "io/io.h"
#include "shape.h"
#define myassert(EX) (void)((EX) || (__assert (#EX, __FILE__, __LINE__),0))

#include <mpi.h>

#define ONROOT int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); if (rank == 0)
int nxexp;
int argnx;
int main() {
  Shape b;
  {
    FILE* q = fopen("nxexp", "r");
    fscanf(q, "%d", &nxexp);
    fclose(q);
  }
  argnx = (1 << nxexp);
  init_grid(argnx);
  FILE* fb = fopen("b.dat", "r");
  b = Read(fb);
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
  boundary({vf});
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
        if (Good(b, x, y, z, Delta))
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
