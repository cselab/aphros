#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>
#include <hdf5_hl.h>
#include <mpi.h>

#include "h5.h"
#include "h5serial.h"

#ifndef H5_HAVE_PARALLEL
#error needs parallel HDF5
#endif

const char* Name = "data";

enum { MAX_SIZE = 4096 };

#define WARN(x)                                     \
  do {                                              \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
    dprint x;                                       \
    fputs("\n", stderr);                            \
  } while (0)
static int dprint(const char* fmt, ...) {
  int r;
  va_list ap;

  va_start(ap, fmt);
  r = vfprintf(stderr, fmt, ap);
  va_end(ap);
  return r;
}

static int h5_strcat(const char a[], const char b[], char* c) {
  while ((*c = *a) != '\0')
    c++, a++;
  while ((*c++ = *b++) != '\0')
    ;
  return 0;
}

static int h5_basename(const char a[], char* b) {
  const char *i, *j;

  for (i = a, j = NULL; *i != '\0'; i++)
    if (*i == '/') j = i;
  if (j == NULL)
    while ((*b++ = *a++) != '\0')
      ;
  else {
    do
      *b++ = *++j;
    while (*j != '\0');
  }
  return 0;
}

static hid_t h5_open(MPI_Comm comm, const char* path) {
  hid_t plist, file;
  char full[MAX_SIZE];

  h5_strcat(path, ".h5", full);
  plist = H5Pcreate(H5P_FILE_ACCESS);
  if (plist < 0) {
    WARN(("H5Pcreate failed"));
    return plist;
  }
  H5Pset_fapl_mpio(plist, comm, MPI_INFO_NULL);
  file = H5Fcreate(full, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
  if (file < 0) {
    WARN(("can't open file '%s'", full));
    return file;
  }
  H5Pclose(plist);
  return file;
}

static herr_t h5_dwrite(
    hid_t file, const char* name, hid_t memspace, hid_t filespace,
    const double* buf) {
  hid_t plist, dset;
  herr_t err;

  dset = H5Dcreate(
      file, name, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT);
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
  err = H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, filespace, plist, buf);
  H5Pclose(plist);
  H5Dclose(dset);
  return err;
}

static int h5_data(
    hid_t file, const int isize[3], const int istart[3], const int iextent[3],
    const double* buf) {
  enum { I, J, K };
  hid_t filespace, memspace;
  herr_t err;
  int dim;

  hsize_t size[] = {isize[I], isize[J], isize[K]};
  hsize_t start[] = {istart[I], istart[J], istart[K]};
  hsize_t extent[] = {iextent[I], iextent[J], iextent[K]};

  dim = 3;
  filespace = H5Screate_simple(dim, size, NULL);
  memspace = H5Screate_simple(dim, extent, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, extent, NULL);

  if (!H5Sselect_valid(filespace)) goto err;
  if (!H5Sselect_valid(memspace)) goto err;
  err = h5_dwrite(file, Name, memspace, filespace, buf);
  H5Sclose(memspace);
  H5Sclose(filespace);
  return err;
err:
  WARN(("invalid selection"));
  WARN(("isize: %d %d %d", isize[I], isize[J], isize[K]));
  WARN(("istart: %d %d %d", istart[I], istart[J], istart[K]));
  WARN(("iextent: %d %d %d", iextent[I], iextent[J], iextent[K]));
  return -1;
}

static int h5_close(hid_t file) {
  return H5Fclose(file);
}

/*
  awk '{gsub(/"/, "\\\""); printf "\"%s\\n\"\\\n", $0}' poc/p.xmf
*/
int h5_xmf(
    const char* path, const char* name, const double origin[3], double spacing,
    const int size[3]) {
  FILE* f;
  enum { I, J, K };
  int i, j, k;
  char base[MAX_SIZE], full[MAX_SIZE], hfile[MAX_SIZE];
  const char s[] =
      "<?xml version=\"1.0\" ?>\n"
      "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
      "<Xdmf Version=\"2.0\">\n"
      " <Domain>\n"
      "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n"
      "     <Topology TopologyType=\"3DCORECTMesh\" Dimensions=\"%d %d %d\"/>\n"
      "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n"
      "       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" "
      "Precision=\"8\" Format=\"XML\">\n"
      "         %.16g %.16g %.16g\n"
      "       </DataItem>\n"
      "       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" "
      "Precision=\"4\" Format=\"XML\">\n"
      "        %.16g %.16g %.16g\n"
      "       </DataItem>\n"
      "     </Geometry>\n"
      "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
      "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" "
      "Precision=\"8\" Format=\"HDF\">\n"
      "        ./%s:/%s\n"
      "       </DataItem>\n"
      "     </Attribute>\n"
      "   </Grid>\n"
      " </Domain>\n"
      "</Xdmf>\n";
  h5_strcat(path, ".xmf", full);
  h5_basename(path, base);
  h5_strcat(base, ".h5", hfile);
  f = fopen(full, "w");
  if (f == NULL) {
    WARN(("can't open file '%s'", full));
    return -1;
  }
  i = size[I];
  j = size[J];
  k = size[K];
  fprintf(
      f, s, i + 1, j + 1, k + 1, origin[I], origin[J], origin[K], spacing,
      spacing, spacing, name, i, j, k, hfile, Name);
  fclose(f);
  return 0;
}

int h5_hdf(
    MPI_Comm comm, const char* path, const int size[3], const int start[3],
    const int extent[3], const double* buf) {
  hid_t file, err;
  int status;

  MPI_Initialized(&status);
  if (status == 0) {
    WARN(("MPI is not initialized"));
    goto err;
  }

  file = h5_open(comm, path);
  if (file < 0) {
    WARN(("h5_open failed for '%s'", path));
    goto err;
  }
  err = h5_data(file, size, start, extent, buf);
  if (err) {
    WARN(("h5_data failed for '%s'", path));
    goto err;
  }
  return h5_close(file);
err:
  MPI_Abort(comm, 2);
  return 1;
}

int h5_silence(void) {
  return H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
}

static void reset_line_end(char* s) {
  for (;;) {
    if (*s == '\0') break;
    if (*s == '\n') {
      *s = '\0';
      break;
    }
    s++;
  }
}

static char* next_line(char* s) {
  for (;; s++) {
    if (*s == '\0') return s;
    if (*s == '\n') return s + 1;
  }
}

static char* to_quote(char* s, char* out) {
  for (;;) {
    if (*s == '\0' || *s == '"') {
      *out = '\0';
      return s;
    }
    *out++ = *s++;
  }
}

int h5_read_xmf(
    const char* path, char* name, double origin[3], double* pspa, int size[3]) {
  enum { X, Y, Z };
  enum { I = Z, J = Y, K = X };

  double spa, spacing[3];
  FILE* f;
  char full[MAX_SIZE], *p, *q, *s;
  size_t n;
  const char* pat;

  h5_strcat(path, ".xmf", full);
  f = fopen(full, "r");
  if (f == NULL) {
    WARN(("can't open '%s'", full));
    goto err;
  }
  if (fseek(f, 0, SEEK_END) != 0) {
    WARN(("fseek failed for '%s'", full));
    goto err;
  }
  n = ftell(f);
  fseek(f, 0, SEEK_SET);
  p = s = malloc(n + 1);
  if (s == NULL) {
    WARN(("malloc(%d) failed for '%s'", n, full));
    goto err;
  }
  if (fread(s, 1, n, f) != n) {
    WARN(("fread failed for '%s'", full));
    goto err;
  }
  s[n++] = '\0';

  p = strstr(p, "<\?xml");
  if (p == NULL) {
    WARN(("not an xml file '%s'", full));
    goto err;
  }

  p = strstr(p, "Xdmf");
  if (p == NULL) {
    WARN(("not an xdmf file '%s'", full));
    goto err;
  }

  p = strstr(p, "\"3DCORECTMesh\"");
  if (p == NULL) {
    WARN(("no 3DCORECTMesh in '%s'", full));
    goto err;
  }

  p = strstr(p, "<DataItem Name=\"Origin\"");
  if (p == NULL) {
    WARN(("no <DataItem Name=\"Origin\" in %s", full));
    goto err;
  }
  q = next_line(p);
  if (sscanf(q, "%lf %lf %lf", &origin[X], &origin[Y], &origin[Z]) != 3) {
    reset_line_end(p);
    WARN(("expecting origin[X] origin[Y] origin[Z] in '%s' after", full));
    WARN(("%s", p));
    goto err;
  }
  p = q;

  p = strstr(p, "<DataItem Name=\"Spacing\"");
  if (p == NULL) {
    WARN(("no <DataItem Name=\"Spacing\""));
    goto err;
  }
  q = next_line(p);
  if (sscanf(q, "%lf %lf %lf", &spacing[X], &spacing[Y], &spacing[Z]) != 3) {
    reset_line_end(p);
    WARN(("expecting spacing[X] spacing[Y] spacing[Z] in '%s' after", full));
    WARN(("%s", p));
    goto err;
  }
  p = q;

  spa = spacing[X];
  if (spacing[Y] != spa || spacing[Z] != spa) {
    WARN(
        ("not unifrom spacing '%.16g %.16g %.16g'", spacing[X], spacing[Y],
         spacing[Z]));
    goto err;
  }

  p = next_line(p);
  pat = "<Attribute Name=\"";
  q = strstr(p, pat);
  q += strlen(pat);
  if (q == NULL) {
    WARN(("no %s", pat));
    goto err;
  }
  q = to_quote(q, name);
  if (*q != '"') {
    reset_line_end(p);
    WARN(("expecting name in '%s' around", full));
    WARN(("%s", p));
    goto err;
  }

  p = next_line(p);
  pat = "<DataItem Dimensions=\"";
  q = strstr(p, pat);
  q += strlen(pat);
  if (q == NULL) {
    WARN(("no %s", pat));
    goto err;
  }
  if (sscanf(q, "%d %d %d", &size[I], &size[J], &size[K]) != 3) {
    reset_line_end(p);
    WARN(("expecting size[I] size[J] size[K] in '%s' around", full));
    WARN(("%s", p));
    goto err;
  }
  free(s);
  fclose(f);
  *pspa = spa;
  return 0;
err:
  return 1;
}

int h5_read_hdf(const char* path, /**/ int size[3], double** pbuf) {
  enum { X, Y, Z };
  enum { I = Z, J = Y, K = X };

  hid_t file;
  herr_t status;
  int rank, n;
  hsize_t dims[99];
  size_t nbytes;
  char full[MAX_SIZE];
  double* buf;

  h5_strcat(path, ".h5", full);
  file = H5Fopen(full, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    WARN(("can't open '%s'", full));
    goto err;
  }
  if (H5LTfind_dataset(file, Name) == 0) {
    WARN(("can't find dataset '%s' in '%s'", Name, full));
    goto err;
  }
  status = H5LTget_dataset_ndims(file, Name, &rank);
  if (status < 0) {
    WARN(("get_dataset_ndims failed for '%s'", full));
    goto err;
  }
  if (rank != 3 && rank != 4) {
    WARN(("rank=%d != (3 or 4) for '%s'", rank, full));
    goto err;
  }
  status = H5LTget_dataset_info(file, Name, dims, NULL, &nbytes);
  if (status < 0) {
    WARN(("get_dataset_info failed for '%s'", full));
    goto err;
  }
  if (nbytes != sizeof(*buf)) {
    WARN(("file: %s", full));
    WARN(("nbytes=%ld != sizeof(*buf)=%d", nbytes, sizeof(*buf)));
    goto err;
  }
  size[X] = dims[I];
  size[Y] = dims[J];
  size[Z] = dims[K];
  n = size[X] * size[Y] * size[Z];
  buf = malloc(n * sizeof(*buf));
  if (buf == NULL) {
    WARN(("can't allocate n = %ld", n));
    goto err;
  }
  status = H5LTread_dataset_double(file, Name, buf);
  if (status < 0) {
    WARN(("read_dataset_double failed for '%s'", full));
    goto err;
  }
  *pbuf = buf;
  return 0;
err:
  return 1;
}

int h5_serial_hdf(const char* path, const int size0[3], const double* buf) {
  enum { X, Y, Z };
  herr_t status;
  hid_t file;
  char full[MAX_SIZE];
  int dimensions = 3;
  const hsize_t size[] = {size0[X], size0[Y], size0[Z]};

  h5_strcat(path, ".h5", full);
  file = H5Fcreate(full, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    WARN(("read_dataset_double failed for '%s'", full));
    goto err;
  }
  status =
      H5LTmake_dataset(file, Name, dimensions, size, H5T_NATIVE_DOUBLE, buf);
  if (status < 0) {
    WARN(("ake_dataset failed for '%s'", full));
    goto err;
  }
  status = H5Fclose(file);
  if (status < 0) {
    WARN(("close failed for '%s'", full));
    goto err;
  }
  return 0;
err:
  return 1;
}
