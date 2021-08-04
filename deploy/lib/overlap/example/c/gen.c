#include <stdio.h>
#include <overlap.h>
#include <gsl/gsl_rng.h>

static double rng(gsl_rng * r, double lo, double hi)
{
  return lo + (hi - lo) * gsl_rng_uniform(r);
}

int main(int args, char **argv) {
  double x, y, z, R, v2, v3, L;
  const gsl_rng_type * T;
  gsl_rng * r;
  long n;
  long i;
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, 123456);
  L = 30;
  n = argv[1] ? atol(argv[1]) : 10;
  for (i = 0; i < n; i++) {
    x = rng(r, -L, L);
    y = rng(r, -L, L);
    z = rng(r, -L, L);
    R = rng(r, 0, L);
    v3 = overlap_3d(x, y, z, R);
    v2 = overlap_2d(x, y, R);
    printf("%.16e %.16e %.16e %.16e %.16e %.16e\n", x, y, z, R, v3, v2);
  }
  gsl_rng_free(r);
}
