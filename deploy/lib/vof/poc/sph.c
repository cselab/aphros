#include <stdio.h>

#include "overlap.h"
#include "vof.h"

int main() {
  // double (*sphere) (double []);
  int i, j, n, dimension;
  double h, a, b, Delta, p[3] = {-0.5, -0.5, -0.5};
  double lo, hi, x, y, z, R;
  R = 10;
  dimension = 3;
  Delta = 1;
  z = 8;
  lo = -10;
  hi = 10;
  n = 10;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      x = lo + (hi - lo) * i / (n - 1);
      y = lo + (hi - lo) * j / (n - 1);
      a = overlap_3d(x, y, z, R);
      b = vof_3d(x, y, z, R);
      printf("%.16e %.16e %.16e %.16e\n", x, y, a, b);
    }
}
