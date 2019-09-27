#include <assert.h>
#include <mpi.h>

#include "fractions.h"
#include ".u/curv/select.h"
#include ".u/io/io.h"
#include ".u/bashape.h"

#define myassert(EX) (void)((EX) || (__assert (#EX, __FILE__, __LINE__),0))

int main() {
  Shape b;
  init_grid(32);
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
}
