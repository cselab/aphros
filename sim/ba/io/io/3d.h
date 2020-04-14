static int io_ver(FILE *f) {
    double a[3], *p;
    foreach_vertex() {
        p = a;
        *p++ = x; *p++ = y; *p++ = z;
        io_write_double(3, a, f);
    }
    return 0;
}

static int io_cell(int nc, scalar m, FILE *f) {
    enum { VOXEL = 8 };
    const int npc = 8;
    const int cell_type = VOXEL;
    int a[npc + 1], *p;
    fprintf(f, "CELLS %d %d\n", nc, (npc + 1) * nc);
    foreach() {
        p = a;
        *p++ = npc;
        *p++ = m[0,0,0]; *p++ = m[1,0,0]; *p++ = m[0,1,0]; *p++ = m[1,1,0];
        *p++ = m[0,0,1]; *p++ = m[1,0,1]; *p++ = m[0,1,1]; *p++ = m[1,1,1];
        io_write_int(sizeof a/sizeof *a, a, f);
    }
    io_write_newline(f);
    fprintf(f, "CELL_TYPES %d\n", nc);
    foreach()
        io_write_int(1, &cell_type, f);
    return 0;
}

