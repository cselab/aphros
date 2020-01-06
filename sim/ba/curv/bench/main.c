#include <assert.h>
#include <mpi.h>

#include "fractions.h"
#include ".u/curv/select.h"
#include ".u/io/io.h"
#include ".u/bashape.h"

#define BDAT "b.dat"

#define myassert(EX) (void)((EX) || (__assert (#EX, __FILE__, __LINE__),0))

int main() {
  Shape b;
  init_grid(16);

  scalar vf[];
  foreach ()  {
    vf[] = 0;
  }

#ifdef BDAT
  const int MAX_SHAPES = 10000;
  int nb = MAX_SHAPES;
  Shape bb[MAX_SHAPES];
  nb = ReadList(BDAT, bb, nb);
  fprintf(stderr, "Reading %d shapes from %s\n", nb, BDAT);
  scalar ft[];
  for (int i = 0; i < nb; ++i) {
    Shape b = bb[i];
    CreateField(b, ft);
    foreach () {
      vf[] = max(vf[], ft[]);
    }
  }
#else
  FILE* fb = fopen("b.dat", "r");
  b = Read(fb);
  fclose(fb);
  origin (0.,0.,0.);
  CreateField(b, vf);
#endif

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

#ifndef NI
  #define NI 10
#endif

  const int ni = NI;
  volatile double a = 0;
  cstats s;
  for (int i = 0; i < ni; ++i) {
    s = curvature(vf, k);
    foreach () {
      a += k[];
      break;
    }
  }
  int nc = s.h + s.f + s.a + s.c; // number of interfacial cells
  printf("cstats: h=%d, f=%d, a=%d, c=%d\n", s.h, s.f, s.a, s.c);
  printf("time elapsed [us/cell]: %.3f\n", kCurvTime / (nc * ni) * 1e6);
}
