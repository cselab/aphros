#include "util.h"
#include "write.h"
#include <mpi.h>
#if dimension == 3
#    include "3d.h"
#else
#    include "2d.h"
#endif

#define assert0( x ) { if( ! (x)) { \
  fprintf(stderr, "%s:%d assertion failed: %s \n", __FILE__, __LINE__, #x);  \
  MPI_Abort(MPI_COMM_WORLD, 1); } }

// Writes fields list for file
// fn: filename
static int iompi(scalar *list, char* fn) {
#   define p(s) fputs(s "\n", f);
    int cr, cs;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &cs);
    MPI_Comm_rank(comm, &cr);

    int root = 0;
    int ir = (cr == root);  // is root

    int nv = 0;  // local vertices
    int nc = 0;  // local cells
    int gnv = 0;  // global vertices
    int gnc = 0;  // global cells

    // local
    foreach_vertex() { nv++; }
    foreach() { nc++; }

    // global
    MPI_Reduce(&nv, &gnv, 1, MPI_INT, MPI_SUM, root, comm);
    MPI_Reduce(&nc, &gnc, 1, MPI_INT, MPI_SUM, root, comm);

    char* bc; // buffer char
    int nb = 0;
    if (ir) {
      nb = max(nb, sizeof(double) * gnv * 3); // points 
      nb = max(nb, sizeof(int) * gnc * 8); // cell types
    } else {
      nb = max(nb, sizeof(double) * nv * 3); // points 
      nb = max(nb, sizeof(int) * nc * 8); // cell types
    }

    bc = malloc(nb);

    FILE* f = NULL;

    if (ir) {
      f = fopen(fn, "w");
      p("# vtk DataFile Version 2.0");
      p("MFER-basilisk");
      p("ASCII");
      p("DATASET UNSTRUCTURED_GRID");
    }

    // points
    {
      double* b = bc;

      ir && fprintf(f, "POINTS %d float\n", gnv);
      int i = 0;
      foreach_vertex() {
        b[i++] = x;
        b[i++] = y;
        b[i++] = z;
      }
      assert0(i == nv * 3);

      if (ir) {
        int rc[cs];  // recvcounts
        int sc = nv * 3;
        MPI_Gather(&sc, 1, MPI_INT, rc, 1, MPI_INT, root, comm);

        int ds[cs];  // displs
        ds[0] = 0;
        for (size_t r = 1; r < cs; ++r) {
          ds[r] = ds[r - 1] + rc[r - 1];
        }
        assert0(ds[cs - 1] + rc[cs - 1] == gnv * 3);

        MPI_Gatherv(MPI_IN_PLACE, nv, MPI_DOUBLE, 
            b, rc, ds, MPI_DOUBLE, root, comm);

        int i = 0;
        for (int v = 0; v < gnv; ++v) {
          fprintf(f, "%g %g %g\n", b[i++], b[i++],  b[i++]);
        }
        assert0(i == gnv * 3);
      } else {
        int sc = nv * 3; // recvcount
        MPI_Gather(&sc, 1, MPI_INT, NULL, 0, MPI_INT, root, comm);
        MPI_Gatherv(b, sc, MPI_DOUBLE, 
            NULL, NULL, NULL, MPI_DOUBLE, root, comm);
      }
    }


    /*
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
    */

    free(bc);

    return 0;
#   undef p
}
