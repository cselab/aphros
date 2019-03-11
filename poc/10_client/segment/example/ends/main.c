#include <stdio.h>
#include <stdlib.h>
#include <segment.h>

#define D (5)

static char **argv0;

int main(int argc, char **argv) {
    enum {AX, AY, BX, BY};
    Scal nx, ny, a;
    Scal ends[4];
    argv0 = argv;
    argv0++;

    nx = atof(*argv0++);
    ny = atof(*argv0++);
    a = atof(*argv0++);

    segment_ends(nx, ny, a, /**/ ends);
    printf("%g %g %g %g\n", ends[AX], ends[AY], ends[BX], ends[BY]);
    return 0;
}
