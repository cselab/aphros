int nx, ny, nz;
double* uu = NULL;

static void ReadField(char* fn) {
  if (!uu) {
    FILE* f = fopen(fn, "r");

    fscanf(f, "%d %d %d", &nx, &ny, &nz);

    uu = malloc(nx * ny * nz * sizeof(double));

    for (int z = 0; z < nz; ++z) {
      for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
          double a;
          fscanf(f, "%lf", &a);
          uu[z * nx * ny + y * nx + x] = a;
        }
      }
    }
    fclose(f);
  }
}

static double ini_t(double x, double y, double z, char* fn) {
  ReadField(fn);

  double h = 1. / nx;
  int ix = MAX(0, MIN(x / h, nx - 1));
  int iy = MAX(0, MIN(y / h, ny - 1));
  int iz = MAX(0, MIN(z / h, nz - 1));
  return uu[iz * nx * ny + iy * nx + ix];
}

static double ini_t2(double x, double y, char* fn) {
  ReadField(fn);

  double h = 1. / nx;
  int ix = MAX(0, MIN(x / h, nx - 1));
  int iy = MAX(0, MIN(y / h, ny - 1));
  return uu[iy * nx + ix];
}

