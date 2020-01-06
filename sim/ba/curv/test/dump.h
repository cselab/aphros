

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
void DumpFacets(scalar c, const char* fn) {
  const int nc = GetNcInter(c); // number of interfacial cells
  const int mm = nc * kMaxFacet;

  coord* xx = (coord*)malloc(mm * sizeof(coord));
  int* pp = (int*)malloc(mm * sizeof(int));
  int* ss = (int*)malloc(mm * sizeof(int));
  int np = 0;
  int nx = 0;

  foreach() {
    if (interfacial (point, c)) {
      coord m = Mycs(point, c);
      double alpha = plane_alpha (c[], m);
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
  WriteVtkPoly(fn, xx, nx, pp, np, ss, "comment", true);

  free(ss);
  free(pp);
  free(xx);
}

// Dumps cross sections of the inteface fragments to vtk.
void DumpLines(scalar c, vector nn, Partstr conf, const char* fn) {
  const int Ns = conf.Ns;

  const int nc = GetNcInter(c); // number of interfacial cells
  const int mm = nc * kMaxSection;

  coord* xx = (coord*)malloc(mm * sizeof(coord));
  int* pp = (int*)malloc(mm * sizeof(int));
  int* ss = (int*)malloc(mm * sizeof(int));
  int nx = 0;
  int np = 0;

  foreach() {
    if (interfacial(point, c)) {
      Trans b = GetPointTrans(point, c, nn);

      for (int s = 0; s < Ns; ++ s) {
        Trans w = GetSectionTrans(s, Ns, b);

        coord ll[kMaxSection];
        int nl = 0;

        Section(point, c, nn, w, ll, &nl);
        GetCrossCurv(point, c, nn, w, conf);

        for (int i = 0; i < nl; ++i) {
          ll[i] = LocToGlb(ll[i], w);
        }

        for (int i = 0; i < nl; i += 2) {
          ss[np] = 2;
          pp[nx] = nx;
          xx[nx] = ll[i];
          ++nx;
          pp[nx] = nx;
          xx[nx] = ll[i + 1];
          ++nx;
          ++np;
        }
      }
    }
  }
  WriteVtkPoly(fn, xx, nx, pp, np, ss, "lines", false);

  free(ss);
  free(pp);
  free(xx);
}
