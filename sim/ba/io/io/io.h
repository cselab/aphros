#include <stdio.h>
#include "geometry.h"
#if dimension == 1
coord mycs (Point point, scalar c) {
  coord n = {1.};
  return n;
}
#elif dimension == 2
# include "myc2d.h"
#else // dimension == 3
# include "myc.h"
#endif
#include "curvature.h"
#include "write.h"
#if dimension == 3
#    include "3d.h"
#else
#    include "2d.h"
#endif

// Set z-copmonent to zero if 2D.
static void A2(coord* p) {
  (void) p;
#if dimension == 2
  p->z = 0;
#endif
}

static coord Add(coord a, coord b) {
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
  return a;
}

static coord Mul(coord a, double k) {
  a.x *= k;
  a.y *= k;
  a.z *= k;
  return a;
}

static coord Mycs(Point point, scalar c) {
  coord m = mycs(point, c);
  A2(&m);
  return m;
}

static int Facets(coord m, double alpha, coord* pp) {
#if dimension == 2
  int nf = facets(m, alpha, pp);
  for (int i = 0; i < nf; ++i) {
    A2(&pp[i]);
  }
  return nf;
#else
  return facets(m, alpha, pp, 1.);
#endif
}

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

const int kMaxFacet = 12; // maximum number of vertices per facet

// Returns the number of interfacial cells.
// vf: volume fraction
static int GetNcInter(scalar vf) {
  int nc = 0;
  foreach () {
    if (interfacial(point, vf)) {
      ++nc;
    }
  }
  return nc;
}

// Writes legacy vtk polydata
// fn: path
// xx: points
// nx: size of xx
// pp: polygons as lists of indices
// np: number of polygons
// ss[i] is size of polygon pp[i]
// cm: comment
// poly: true: polygons, false: lines
void WriteVtkPoly(const char* fn, coord* xx, int nx,
                  int* pp, int np, int* ss, const char* cm, bool poly) {
  FILE* o = fopen(fn, "w");
  fprintf(o, "# vtk DataFile Version 2.0\n");
  fprintf(o, "%s\n", cm);

  fprintf(o, "ASCII\n");
  fprintf(o, "DATASET POLYDATA\n");

  fprintf(o, "POINTS %d float\n", nx);
  for (int i = 0; i < nx; ++i) {
    coord x = xx[i];
    fprintf(o, "%g %g %g\n", x.x, x.y, x.z);
  }

  int na = 0; // total number of vortices
  for (int i = 0; i < np; ++i) {
    na += ss[i];
  }
  fprintf(o, "%s %d %d\n", poly ? "POLYGONS" : "LINES", np, np + na);
  int k = 0;
  for (int i = 0; i < np; ++i) {
    fprintf(o, "%d", ss[i]);
    for (int j = 0; j < ss[i]; ++j) {
      fprintf(o, " %d", pp[k++]);
    }
    fprintf(o, "\n");
  }

  fclose(o);
}

// Dumps interface fragments to vtk.
// vf: volume fraction
void DumpFacets(scalar vf, const char* filename) {
  const int nc = GetNcInter(vf); // number of interfacial cells
  const int mm = nc * kMaxFacet;

  coord* xx = (coord*)malloc(mm * sizeof(coord));
  int* pp = (int*)malloc(mm * sizeof(int));
  int* ss = (int*)malloc(mm * sizeof(int));
  int np = 0;
  int nx = 0;

  foreach() {
    if (interfacial (point, vf)) {
      coord m = Mycs(point, vf);
      double alpha = plane_alpha (vf[], m);
      coord p[kMaxFacet];
      int nf = Facets(m, alpha, p);

      coord r = {x, y, z};

      ss[np] = nf;
      for (int i = 0; i < nf; ++i) {
        pp[nx] = nx;
        xx[nx] = Add(r, Mul(p[i], Delta));
        ++nx;
      }
      ++np;
    }
  }
  WriteVtkPoly(filename, xx, nx, pp, np, ss, "comment", true);

  free(ss);
  free(pp);
  free(xx);
}
