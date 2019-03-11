#include <stdio.h>
#include <stdlib.h>
#include <segment.h>

#define D (5)

static char **argv0;

Scal tof(char **s) {
    if (s == NULL) {
        fprintf(stderr, "not enough argument\n");
        exit(2);
    }
    return atof(*s);
}

int main(int argc, char **argv) {
    enum {AX, AY, BX, BY};
    Scal nx, ny, a;
    Scal ends[4];

    argv++;
    nx = tof(argv++);
    ny = tof(argv++);
    a = tof(argv++);

    segment_ends(nx, ny, a, /**/ ends);
    printf("%g %g %g %g\n", ends[AX], ends[AY], ends[BX], ends[BY]);
    return 0;
}
