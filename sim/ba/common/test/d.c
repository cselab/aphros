#include <tgmath.h>

#include "div.h"
#include <stdio.h>

double A(
    double u0, double u1, double u2,
    double u3, double u4, double u5,
    double u6, double u7, double u8
    ) {
  double u[DIM] = {u0, u1, u2, u3, u4, u5, u6, u7, u8};
  return Curv0(u);
}

int main() {
  printf("%g , %g\n", A(0, 0, 0, 0.5, 0.5, 0.5, 1, 1, 1), 0.);
  printf("%g , %g\n", A(0, 0, 0, 0.2, 0.2, 0.2, 1, 1, 1), 0.);
}
