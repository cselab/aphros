static int io_write_int(int n, const int *a, FILE *f) {
    int i = 0;
    fprintf(f, "%d", a[i++]);
    while (i < n)
        fprintf(f, "  %d", a[i++]);
    fputs("\n", f);
    return 0;
}

static int io_write_double(int n, const double *a, FILE *f) {
    int i;
    for (i = 0; i < n; i++) {
        if (i > 0)
            fputs(" ", f);
        fprintf(f, "%.16e", a[i]);
    }
    fputs("\n", f);
    return 0;
}

static int io_write_type(FILE *f) { fputs("ASCII\n", f); }
static int io_write_newline(FILE *f) { return 0; }
