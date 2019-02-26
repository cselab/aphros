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
    {
      int nv_ = (ir ? gnv : nv);
      int nc_ = (ir ? gnc : nc);
      nb = max(nb, sizeof(double) * nv_ * 3); // points 
      nb = max(nb, sizeof(int) * nc_ * 9); // cell nodes
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

      int i = 0;
      foreach_vertex() {
        b[i++] = z;
        b[i++] = y;
        b[i++] = x;
      }
      assert0(i == nv * 3);

      if (ir) {
        fprintf(f, "POINTS %d float\n", gnv);
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

        fprintf(f, "\n");
      } else {
        int sc = nv * 3; // recvcount
        MPI_Gather(&sc, 1, MPI_INT, NULL, 0, MPI_INT, root, comm);
        MPI_Gatherv(b, sc, MPI_DOUBLE, 
            NULL, NULL, NULL, MPI_DOUBLE, root, comm);
      }
    }
    
    // XXX: adhoc 3D

    // cell nodes
    {
      int* b = bc;

      const int npc = 8; // points per cell
      const int ic = 9; // offset per cell

      vertex scalar m[];
      {
        int i = 0;
        foreach_vertex() { m[] = i++; }
      }

      int i = 0;
      foreach() {
        b[i++] = npc;
        b[i++] = m[0,0,0]; b[i++] = m[1,0,0]; 
        b[i++] = m[0,1,0]; b[i++] = m[1,1,0];
        b[i++] = m[0,0,1]; b[i++] = m[1,0,1]; 
        b[i++] = m[0,1,1]; b[i++] = m[1,1,1];
      }
      assert0(i == nc * ic);

      if (ir) {
        fprintf(f, "CELLS %d %d\n", gnc, gnc * ic);

        int rc[cs];  // recvcounts
        int sc = nc * ic;
        MPI_Gather(&sc, 1, MPI_INT, rc, 1, MPI_INT, root, comm);

        int nnv[cs];  // points on ranks
        MPI_Gather(&nv, 1, MPI_INT, nnv, 1, MPI_INT, root, comm);


        int ds[cs];  // displs
        ds[0] = 0;
        for (size_t r = 1; r < cs; ++r) {
          ds[r] = ds[r - 1] + rc[r - 1];
        }
        assert0(ds[cs - 1] + rc[cs - 1] == gnc * ic);

        MPI_Gatherv(MPI_IN_PLACE, nc, MPI_INT, 
            b, rc, ds, MPI_INT, root, comm);

        int i = 0;
        int cb = 0; // base cell
        for (int r = 0; r < cs; ++r) {
          for (int c = 0; c < rc[r] / ic; ++c) {
            fprintf(f, "%d", b[i++]);
            for (int j = 1; j < ic; ++j) {
              fprintf(f, " %d", b[i++] + cb);
            }
            fprintf(f, "\n");
          }
          cb += nnv[r];
        }
        assert0(i == gnc * ic);
        fprintf(f, "\n");

        const int ct = 11; // cell type, voxel
        fprintf(f, "CELL_TYPES %d\n", gnc);
        for (int c = 0; c < gnc; ++c) {
          fprintf(f, "%d \n", ct);
        }
        fprintf(f, "\n");
      } else {
        int sc = nc * ic; // recvcount
        MPI_Gather(&sc, 1, MPI_INT, NULL, 0, MPI_INT, root, comm);
        MPI_Gather(&nv, 1, MPI_INT, NULL, 0, MPI_INT, root, comm);
        MPI_Gatherv(b, sc, MPI_INT, 
            NULL, NULL, NULL, MPI_INT, root, comm);
      }
    }

    ir && fprintf(f, "CELL_DATA %d\n", gnc);

    for (scalar s in list) {
      double* b = bc;

      int i = 0;
      foreach() {
        b[i++] = s[];
      }
      assert0(i == nc);

      if (ir) {
        fprintf(f, "SCALARS %s float\n", s.name);
        p("LOOKUP_TABLE default");

        int rc[cs];  // recvcounts
        int sc = nc;
        MPI_Gather(&sc, 1, MPI_INT, rc, 1, MPI_INT, root, comm);

        int ds[cs];  // displs
        ds[0] = 0;
        for (size_t r = 1; r < cs; ++r) {
          ds[r] = ds[r - 1] + rc[r - 1];
        }
        assert0(ds[cs - 1] + rc[cs - 1] == gnc);

        MPI_Gatherv(MPI_IN_PLACE, nc, MPI_DOUBLE, 
            b, rc, ds, MPI_DOUBLE, root, comm);

        int i = 0;
        for (int c = 0; c < gnc; ++c) {
          fprintf(f, "%.10g ", b[i++]);
          fprintf(f, "\n");
        }
        assert0(i == gnc);

        fprintf(f, "\n");
      } else {
        int sc = nc; // recvcount
        MPI_Gather(&sc, 1, MPI_INT, NULL, 0, MPI_INT, root, comm);
        MPI_Gatherv(b, sc, MPI_DOUBLE, 
            NULL, NULL, NULL, MPI_DOUBLE, root, comm);
      }
    }

    free(bc);

    ir && fclose(f);

    return 0;
#   undef p
}
