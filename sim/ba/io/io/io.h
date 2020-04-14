#include <stdio.h>
#include "write.h"
#if dimension == 3
#    include "3d.h"
#else
#    include "2d.h"
#endif

static int io_write(scalar *list, FILE *f) {
    int nv, nc;
    vertex scalar m[];
    nv = nc = 0;
    foreach_vertex()
        m[] = nv++;
    foreach()
        nc++;
    if (fputs("# vtk DataFile Version 2.0\n", f) == EOF) {
        fprintf(stderr, "io: fputs failed\n");
        return 1;
    }
    fputs("aphros-basilisk\n", f);
    io_write_type(f);
    fputs("DATASET UNSTRUCTURED_GRID\n", f);
    fprintf(f, "POINTS %d float\n", nv);
    io_ver(f);
    io_write_newline(f);
    io_cell(nc, m, f);
    io_write_newline(f);
    fprintf(f, "CELL_DATA %d\n", nc);
    for (scalar s in list) {
        fprintf(f, "SCALARS %s float\n", s.name);
        fputs("LOOKUP_TABLE default\n", f);
        foreach()
            io_write_double(1, &s[], f);
        io_write_newline(f);
    }
    return 0;
}

