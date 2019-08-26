#include <h5.h>

#define BAH5_WARN(x) do {                                       \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__);             \
    bah5_dprint x;                                              \
    fputs("\n", stderr);                                        \
  } while (0)

static int bah5_dprint(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
}

static int
bah5_write(scalar s, const char *name, const char *path)
{
  enum {X, Y, Z};
  enum {I = Z, J = Y, K = X};
  int status;
  int e, nbuf, ext[3], siz[3], sta[3];
  double spa, ori[3], *buf;
  int i, j, k, idx, nc;

  if (!grid) {
    BAH5_WARN(("grid is not initialized"));
    goto err;
  }
  if (depth() != grid->maxdepth)
    BAH5_WARN(("not uniform grid"));
  MPI_Initialized(&status);
  if (!status)
    BAH5_WARN(("MPI is not initialized"));

  e = 1 << depth();
  ext[I] = ext[J] = ext[K] = e;
  siz[I] = e*mpi_dims[X];
  siz[J] = e*mpi_dims[Y];
  siz[K] = e*mpi_dims[Z];

  sta[I] = e*mpi_coords[X];
  sta[J] = e*mpi_coords[Y];
  sta[K] = e*mpi_coords[Z];

  ori[X] = X0; ori[Y] = Y0; ori[Z] = Z0;
  if (siz[X] == 0) {
    BAH5_WARN(("siz[X] == 0"));
    goto err;
  }
  spa = L0/siz[X];
  nbuf = e*e*e;
  buf = malloc(nbuf*sizeof(*buf));
  if (buf == NULL) {
    BAH5_WARN(("fail to allocate 'nbuf = %d'", nbuf));
    goto err;
  }
  nc = 0;
  foreach() {
    i = point.i - GHOSTS;
    j = point.j - GHOSTS;
    k = point.k - GHOSTS;
    idx = k*e*e + j*e + i;
    buf[idx] = s[];
    if (idx >= nbuf) {
      BAH5_WARN(("idx=%d >= nbuf=%d", idx, nbuf));
      goto err;
    }
    nc++;
  }
  if (nc != nbuf) {
    BAH5_WARN(("nc=%d  != nbuf=%d", nc, nbuf));
    goto err;
  }
  status = h5_xmf(path, name, ori, spa, siz);
  if (status != 0) {
    BAH5_WARN(("can't write xmf '%s'", path));
    goto err;
  }
  status = h5_hdf(MPI_COMM_WORLD, path, siz, sta, ext, buf);
  if (status != 0) {
    BAH5_WARN(("can't write hdf '%s'", path));
    goto err;
  }
  free(buf);
  return 0;
err:
  MPI_Abort(MPI_COMM_WORLD, 1);
}

static int
bah5_list(scalar *list, const char *templ)
{
  char path[4096];
  for (scalar s in list) {
    if (snprintf(path, 4096, templ, s.name) < 0) {
      BAH5_WARN(("no '%%s' in templ = '%s'", templ));
      return 0;
    }
    bah5_write(s, s.name, path);
  }
}

static int
bah5_list0(scalar *list, char *name[], const char *templ)
{
  char path[4096];
  for (scalar s in list) {
    if (snprintf(path, 4096, templ, *name) < 0) {
      BAH5_WARN(("no '%%s' in templ = '%s'", templ));
      return 0;
    }
    bah5_write(s, *name, path);
    name++;
  }
}
