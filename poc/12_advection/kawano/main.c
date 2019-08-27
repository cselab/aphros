#include <stdio.h>
#include "sz.h"

int main() {
    double alpha, alpha0, a, b, c, v, g, c0;
    alpha = 0.8;
    a = b = c = 1/3.0;
    g = 0.2;
    c0 = 0.35;
    v = calc_v(alpha, a, b, c);
    alpha0 = calc_alpha(v, a, b, c);
    printf("v = %g\n", v);
    printf("%g %g %g\n", alpha0, alpha, alpha0 - alpha);
    printf("%g\n", calc_flux_plic(g, c0, a, b, c));
}
