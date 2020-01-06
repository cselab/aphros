#include <assert.h>
#include <mpi.h>

#include "fractions.h"
#include ".u/curv/select.h"
#include ".u/io/io.h"
#include ".u/bashape.h"
#include ".u/bah5.h"

#define myassert(EX) (void)((EX) || (__assert (#EX, __FILE__, __LINE__),0))

double levelset(double x, double y, double z) {
  double r = sqrt(sq(x - 0.5) + sq(y - 0.5) + sq(z - 0.5));
  r *= 1 + 0.3 * sin(x * 10) * sin(y * 15);
  return r - 0.3;
}

int main() {
  Shape b;
  init_grid(32);
  origin (0.,0.,0.);
  scalar vf[];
  fraction(vf, levelset(x, y, z));

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
    FILE* q = fopen("curv", "w");
    foreach() {
      if (vf[] > 0. && vf[] < 1.) {
        if (Good(b, x, y, z, Delta))
          fprintf(q, "%.16g\n", k[] * kc);
      }
    }
    fclose(q);
  }

  {
    bah5_list({vf}, "%s");
  }
}
