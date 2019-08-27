#include <stdio.h>
double calc_v(double alpha, double, double, double);
double calc_alpha(double v, double, double, double);

int main() {
    double alpha, alpha0, a, b, c, v;
    alpha = 0.8;
    a = b = c = 1/3.0;
    v = calc_v(alpha, a, b, c);
    alpha0 = calc_alpha(v, a, b, c);
    printf("v = %g\n", v);
    printf("%g %g %g\n", alpha0, alpha, alpha0 - alpha);
}
