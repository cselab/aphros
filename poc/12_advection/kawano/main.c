#include <stdio.h>
double calc_v(double alpha, double, double, double);

int main() {
    double alpha, a, b, c, v;
    
    alpha = 0.5;
    a = b = c = 1/3.0;

    v = calc_v(alpha, a, b, c);
    printf("%g\n", v);
}
