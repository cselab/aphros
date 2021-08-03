#include <stdio.h>
#include <overlap.h>
#include <gsl/gsl_rng.h>

static double rng(gsl_rng * r, double lo, double hi)
{
  return (hi - lo) * gsl_rng_uniform(r) - lo;
}

int main(int args, char **argv) {
  double x, y, z, R, v, L;
  const gsl_rng_type * T;
  gsl_rng * r;
  int n;
  int i;
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, 123456);
  L = 10;
  n = argv[1] ? atoi(argv[1]) : 10;
  for (i = 0; i < n; i++) {
    x = rng(r, -L, L);
    y = rng(r, -L, L);
    z = rng(r, -L, L);
    R = rng(r, -L, L);
    v = overlap_3d(x, y, z, R);
    printf("%.16e %.16e %.16e %.16e %.16e\n", x, y, z, R, v);
  }
  gsl_rng_free(r);
}
