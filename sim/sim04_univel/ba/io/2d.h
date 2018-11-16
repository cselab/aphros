static int vtk_ver(FILE *f) {
    double a[3], *p;
    foreach_vertex() {
        p = a;
        *p++ = x; *p++ = y; *p++ = 0;
        write_double(3, a, f);
    }
    return 0;
}

static int vtk_cell(int nc, scalar m, FILE *f) {
    const int npc = 4;
    const int cell_type = PIXEL;
    int a[npc + 1], *p;
    fprintf(f, "CELLS %d %d\n", nc, (npc + 1) * nc);
    foreach() {
        p = a;
        *p++ = npc;
        *p++ = m[0,0]; *p++ = m[1,0]; *p++ = m[0,1]; *p++ = m[1,1];
        write_int(SIZE(a), a, f);
    }
    write_newline(f);
    fprintf(f, "CELL_TYPES %d\n", nc);
    foreach()
        write_int(1, &cell_type, f);
    return 0;
}

