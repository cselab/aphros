#include <stdio.h>
#include <segment.h>

#define D (5)

int main() {
    int i, j;
    double alpha[D*D], *normal, *a;
    for (i = 0; i < D; i++)
        for (j = 0; j < D; j++)
            alpha[D*i + j] = (i*j + 0.0)/(D*D);

    segment_get(alpha, &normal, &a);

    for (i = 0; i < D*D; i++)
        printf("%g %g\n", normal[i], normal[i + 2*D]);

    for (i = 0; i < D*D; i++)
        printf("%g\n", a[i]);
}
