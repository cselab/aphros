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
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
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
h5_dwrite(hid_t file, const char *name, hid_t memspace, hid_t filespace, double *buf)
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
h5_data(hid_t file, int *isize, int *istart, int *iextent, double *buf)
{
  hid_t filespace, memspace;
  herr_t err;
  int dim;
  hsize_t size[] = {
    isize[0], isize[1], isize[2]};
  hsize_t start[] = {
    istart[0], istart[1], istart[2]};
  hsize_t extent[] = {
    iextent[0], iextent[1], iextent[2]};

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
  WARN(("isize: %d %d %d", isize[0], isize[1], isize[2]));
  WARN(("istart: %d %d %d", istart[0], istart[1], istart[2]));
  WARN(("iextent: %d %d %d", iextent[0], iextent[1], iextent[2]));
  abort();
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
h5_xmf(const char *path, const char *name, double origin[3], double spacing, int size[3])
{
  FILE *f;
  enum {X, Y, Z};
  int x, y, z;
  char base[MAX_SIZE], full[MAX_SIZE], hfile[MAX_SIZE];
  const char s[] =
"<?xml version=\"1.0\" ?>\n"\
"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"\
"<Xdmf Version=\"2.0\">\n"\
" <Domain>\n"\
"   <Grid GridType=\"Uniform\">\n"\
"     <Time Value=\"0.000000e+00\"/>\n"\
"     <Topology TopologyType=\"3DCORECTMesh\" Dimensions=\"%d %d %d\"/>\n"\
"     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n"\
"       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n"\
"         %.16g %.16g %.16g\n"\
"       </DataItem>\n"\
"       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n"\
"        %.16g %.16g %.16g\n"\
"       </DataItem>\n"\
"     </Geometry>\n"\
"     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n"\
"       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"\
"        %s:/data\n"\
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
  x = size[X]; y = size[Y]; z = size[Z];
  fprintf(f, s, x, y, z, origin[X], origin[Y], origin[Z], spacing, spacing, spacing, name, x, y, z, hfile);
  fclose(f);
  return 0;
}

int h5_hdf(MPI_Comm comm, const char *path, int size[3], int start[3], int extent[3], double *buf)
{
  hid_t file;
  file = h5_open(comm, path);
  if (file < 0) {
    WARN(("h5_open failed for '%s'", path));
    return 1;
  }
  h5_data(file, size, start, extent, buf);
  return h5_close(file);
}

int
h5_silence(void)
{
  return H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
}
