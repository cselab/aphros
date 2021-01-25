/* c99 eig_values.c -I.. -L.. -lchmath `pkg-config --libs gsl` */
#include <stdio.h>
#include <chmath.h>

int main() {
  double A[6] = {1, 2, 3, 4, 5, 6};
  double w[3];
  chmath_eig_values(A, w);
  printf("%g %g %g\n", w[0], w[1], w[2]);
}
