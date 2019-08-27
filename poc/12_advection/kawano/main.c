#include <stdio.h>
#include "sz.h"

int main() {
    double al, al0, a, b, c, v, g, c0;
    al = 0.8;
    a = b = c = 1/3.0;
    g = 0.2;
    c0 = 0.35;
    v = vof(al, a, b, c);
    al0 = alpha(v, a, b, c);
    printf("v = %g\n", v);
    printf("%g %g %g\n", al, al0, al0 - al);
    printf("%g\n", flux_plic(g, c0, a, b, c));
}
