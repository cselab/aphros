#include <stdio.h>
#include <segment.h>

#define D (5)

int main() {
    int i, j;
    double alpha[D*D], *normal, *a;
    for (j = 0; j < D; j++)
        for (i = 0; i < D; i++)
            alpha[D*j + i] = (i + 0.0)/(D + 1);

    segment_get(alpha, &normal, &a);

    for (i = 0; i < D*D; i++)
        printf("%g %g\n", normal[i], normal[i + D*D]);

//    for (i = 0; i < D*D; i++)
//        printf("%g\n", a[i]);
}
