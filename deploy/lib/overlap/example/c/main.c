#include <stdio.h>

#include <overlap.h>

int main() {
  double x, y, z, r, v;

  x = y = z = 0;
  r = 0.1;
  v = overlap_3d(x, y, z, r);
  printf("%g\n", v);
}
