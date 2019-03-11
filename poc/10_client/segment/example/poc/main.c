#include <stdio.h>
#include <tgmath.h>
#include <segment.h>

#define D (5)

int main(void) {
#define sq(x) ((x)*(x))
    int i, j;
    double alpha[D*D], *normal, *a, *s, nx, ny;

    for (j = 0; j < D; j++)
        for (i = 0; i < D; i++) {
            alpha[D*j + i] = 1/(sq(i-2.0) + sq(j-2.0) + 1);
        }

    segment_get(alpha, &normal, &a, &s);

    for (i = 0; i < D*D; i++) printf("%g %g\n", normal[2*i], normal[2*i + 1]);
    for (i = 0; i < D*D; i++) printf("%g\n", a[i]);
    for (i = 0; i < D*D; i++) printf("%g %g %g %g\n", s[4*i], s[4*i + 1], s[4*i + 2], s[4*i + 3]);

    for (i = 0; i < D*D; i++) {
        printf("%g %g\n", s[4*i], s[4*i + 1]);
        printf("%g %g\n", s[4*i + 2], s[4*i + 3]);
        printf("\n");
    }

    segment_norm(2, 2, alpha, &nx, &ny);

    i = 2*D + 2;
    //printf("%.16e %.16e\n", normal[2*i], normal[2*i + 1]);
    //printf("%.16e %.16e\n", nx, ny);
}
