#include <stdio.h>
#include <tgmath.h>
#include <segment.h>

#define D (5)

int main() {
    int i, j;
    double alpha[D*D], *normal, *a, nx, ny;

    for (j = 0; j < D; j++)
        for (i = 0; i < D; i++)
            alpha[D*j + i] = (cos(i)*sin(j*i) + 0.0)/(D + 1);

    segment_get(alpha, &normal, &a);

    for (i = 0; i < D*D; i++)
        printf("%g %g\n", normal[2*i], normal[2*i + 1]);

    for (i = 0; i < D*D; i++)
        printf("%g\n", a[i]);

    segment_norm(2, 2, alpha, &nx, &ny);

    i = 2*D + 2;
    printf("%.16e %.16e\n", normal[2*i], normal[2*i + 1]);
    printf("%.16e %.16e\n", nx, ny);
}
