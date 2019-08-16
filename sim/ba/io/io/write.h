#include <stdio.h>

static int write_int(int n, const int *a, FILE *f) {
    int i = 0;
    fprintf(f, "%d", a[i++]);
    while (i < n)
        fprintf(f, "  %d", a[i++]);
    fputs("\n", f);
    return 0;
}

static int write_double(int n, const double *a, FILE *f) {
    int i = 0;
    fprintf(f, "%.16g", a[i++]);
    while (i < n)
        fprintf(f, "  %.16g", a[i++]);
    fputs("\n", f);
    return 0;
}

static int write_type(FILE *f) { fputs("ASCII\n", f); }
static int write_newline(FILE *f) { return 0; }
