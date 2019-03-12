#include <stdio.h>
#include <stdlib.h>
#include <segment.h>

#define D (5)

static char **argv0;

Scal tof(char *s) {
    if (s == NULL) {
        fprintf(stderr, "not enough argument\n");
        exit(2);
    }
    return atof(s);
}

int main(int argc, char **argv) {
    Scal nx, ny, u;
    Scal ends[4];

    argv++;
    nx = tof(*argv++);
    ny = tof(*argv++);
    u = tof(*argv++);

    printf("%g\n", segment_line(nx, ny, u));
    return 0;
}
