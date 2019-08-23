#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <hdf5.h>

#include "h5.h"

#ifndef H5_HAVE_PARALLEL
#error needs parallel HDF5
#endif

#define MAX_SIZE 4096
#define WARN(x) do {                                    \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__);     \
    dprint x;                                           \
    fputs("\n", stderr);                                \
  } while (0)
static int dprint(const char *fmt, ...) {
  int r;
  va_list ap;
  va_start(ap, fmt);
  r = vfprintf(stderr, fmt, ap);
  va_end(ap);
  return r;
}

static int
h5_strcat(const char a[], const char b[], char *c)
{
  while ((*c = *a) != '\0') c++, a++;
  while ((*c++ = *b++) != '\0');
  return 0;
}

static int
h5_basename(const char a[], char *b)
{
  const char *i, *j;
  for (i = a, j = NULL; *i != '\0'; i++)
    if (*i == '/') j = i;
  if (j == NULL)
    while ((*b++ = *a++) != '\0');
  else {
    do *b++ = *++j;
    while (*j != '\0');
  }
  return 0;
}

static hid_t
h5_open(MPI_Comm comm, const char *path)
{
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

static herr_t
h5_dwrite(hid_t file, const char *name, hid_t memspace, hid_t filespace, const double *buf)
{
  hid_t plist, dset;
  herr_t err;

  dset = H5Dcreate(file, name, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
  err = H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, filespace, plist, buf);
  H5Pclose(plist);
  H5Dclose(dset);
  return err;
}

static int
h5_data(hid_t file, const int isize[3], const int istart[3], const int iextent[3], const double *buf)
{
  enum {I, J, K};
  hid_t filespace, memspace;
  herr_t err;
  int dim;
  hsize_t size[] = {
    isize[I], isize[J], isize[K]};
  hsize_t start[] = {
    istart[I], istart[J], istart[K]};
  hsize_t extent[] = {
    iextent[I], iextent[J], iextent[K]};

  dim = 3;
  filespace = H5Screate_simple(dim, size, NULL);
  memspace = H5Screate_simple(dim, extent, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, extent, NULL);

  if (!H5Sselect_valid(filespace)) goto err;
  if (!H5Sselect_valid(memspace)) goto err;
  err = h5_dwrite(file, "data", memspace, filespace, buf);
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

static int
h5_close(hid_t file)
{
  return H5Fclose(file);
}

/*
  awk '{gsub(/"/, "\\\""); printf "\"%s\\n\"\\\n", $0}' poc/p.xmf
*/
int
h5_xmf(const char *path, const char *name, const double origin[3], double spacing, const int size[3])
{
  FILE *f;
  enum {I, J, K};
  int i, j, k;
  char base[MAX_SIZE], full[MAX_SIZE], hfile[MAX_SIZE];
  const char s[] =
"<?xml version=\"1.0\" ?>\n"\
"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"\
"<Xdmf Version=\"2.0\">\n"\
" <Domain>\n"\
"   <Grid GridType=\"Uniform\">\n"\
"     <Topology TopologyType=\"3DCORECTMesh\" Dimensions=\"%d %d %d\"/>\n"\
"     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n"\
"       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n"\
"         %.16g %.16g %.16g\n"\
"       </DataItem>\n"\
"       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n"\
"        %.16g %.16g %.16g\n"\
"       </DataItem>\n"\
"     </Geometry>\n"\
"     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n"\
"       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"\
"        ./%s:/data\n"\
"       </DataItem>\n"\
"     </Attribute>\n"\
"   </Grid>\n"\
" </Domain>\n"\
"</Xdmf>\n"\
    ;
  h5_strcat(path, ".xmf", full);
  h5_basename(path, base);
  h5_strcat(base, ".h5", hfile);
  f = fopen(full, "w");
  if (f == NULL) {
    WARN(("can't open file '%s'", full));
    return -1;
  }
  i = size[I]; j = size[J]; k = size[K];
  fprintf(f, s, i + 1, j + 1, k + 1, origin[I], origin[J], origin[K], spacing, spacing, spacing, name, i, j, k, hfile);
  fclose(f);
  return 0;
}

int h5_hdf(MPI_Comm comm, const char *path, const int size[3], const int start[3], const int extent[3], const double *buf)
{
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

int
h5_silence(void)
{
  return H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
}
