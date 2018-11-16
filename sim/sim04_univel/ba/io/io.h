#include "util.h"
#include "write.h"
#if dimension == 3
#    include "3d.h"
#else
#    include "2d.h"
#endif

static int io(scalar *list, FILE *f) {
#   define p(s) fputs(s "\n", f);
    int nv, nc;
    vertex scalar m[];
    nv = nc = 0;
    foreach_vertex() m[] = nv++;
    foreach() nc++;
    p("# vtk DataFile Version 2.0");
    p("MFER-basilisk");
    write_type(f);
    p("DATASET UNSTRUCTURED_GRID");
    fprintf(f, "POINTS %d float\n", nv);
    vtk_ver(f);
    write_newline(f);
    vtk_cell(nc, m, f);
    write_newline(f);
    fprintf(f, "CELL_DATA %d\n", nc);
    for (scalar s in list) {
        fprintf(f, "SCALARS %s float\n", s.name);
        p("LOOKUP_TABLE default");
        foreach()
            write_double(1, &s[], f);
        write_newline(f);
    }
    return 0;
#   undef p
}
