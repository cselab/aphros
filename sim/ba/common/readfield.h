int ReadField(scalar c, char* fn) {
  FILE* f = fopen(fn, "r");

  if (!f) {
    return 0;
  }

  int nx, ny, nz;
  fscanf(f, "%d %d %d", &nx, &ny, &nz);

  myassert(nx == argnx);
  myassert(ny == argnx);

#if dimension == 2
  myassert(nz == 1);
#elif dimension == 3
  myassert(nz == argnx);
#endif

  double* uu = malloc(sizeof(double) * nx * ny * nz);

  for (int z = 0; z < nz; ++z) {
    for (int y = 0; y < ny; ++y) {
      for (int x = 0; x < nx; ++x) {
        double a;
        fscanf(f, "%lf", &a);
        uu[z * ny * nx + y * nx + x] = a;
      }
    }
  }
  fclose(f);

  int i = 0;
  foreach() {
    double h = Delta;
    double hmin = 1. / argnx;
    int ix = max(0, min(x / hmin, nx - 1));
    int iy = max(0, min(y / hmin, ny - 1));
    int iz = max(0, min(z / hmin, nz - 1));
    c[] = uu[iz * ny * nx + iy * nx + ix];
  }

  boundary ({c});

  free(uu);

  return 1;
}
