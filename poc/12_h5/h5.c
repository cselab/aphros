#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <hdf5.h>

#include "h5.h"

#ifndef H5_HAVE_PARALLEL
#error needs parallel HDF5
#endif

#define MAX_SIZE 4096

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

hid_t
h5_open(MPI_Comm comm, const char *path)
{
	hid_t plist, file;
	char full[MAX_SIZE];

	h5_strcat(path, ".h5", full);
	plist = H5Pcreate(H5P_FILE_ACCESS);
	if (plist < 0) {
		fprintf(stderr, "%s:%d: H5Pcreate failed\n", __FILE__, __LINE__);
		return plist;
	}
	H5Pset_fapl_mpio(plist, comm, MPI_INFO_NULL);
	file = H5Fcreate(full, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
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

int
h5_data(hid_t file, int *isize, int *istart, int *iextent, double *buf)
{
	hid_t filespace, memspace;
	herr_t err;
	int dim;
	hsize_t size[] = {
		isize[0], isize[1], isize[2], 1	};
	hsize_t start[] = {
		istart[0], istart[1], istart[2], 0	};
	hsize_t extent[] = {
		iextent[0], iextent[1], iextent[2], 0	};

	dim = 4;
	filespace = H5Screate_simple(dim, size, NULL);
	memspace = H5Screate_simple(dim, extent, NULL);
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, extent, NULL);
	err = h5_dwrite(file, "data", memspace, filespace, buf);
	H5Sclose(memspace);
	H5Sclose(filespace);
	return err;
}

static int
h5_swrite(hid_t file, int root, const char *name, int isize, double h)
{
	hid_t space, dset;
	herr_t err;
	int i, dim;
	double *buf;
	hsize_t size[] = {
		isize + 1		};
	buf = malloc(sizeof(buf[0])*(isize + 1));
	if (buf == NULL) {
		fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
		return -1;
	}
	for (i = 0; i < isize + 1; i++)
		buf[i] = i * h;
	dim = 1;
	space = H5Screate_simple(dim, size, NULL);
	if (space < 0) {
		fprintf(stderr, "%s:%d: H5Screate_simple\n", __FILE__, __LINE__);
		return -1;
	}
	dset = H5Dcreate(file, name, H5T_NATIVE_DOUBLE,	 space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dset < 0) {
		fprintf(stderr, "%s:%d: H5Dcreate failed\n", __FILE__, __LINE__);
		return -1;
	}

	if (!root) H5Sselect_none(space);
	err = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
	if (err < 0) {
		fprintf(stderr, "%s:%d: H5Dwrite failed\n", __FILE__, __LINE__);
		return -1;
	}
	H5Sclose(space);
	H5Dclose(dset);
	free(buf);
	return 0;
}

int
h5_spacing(hid_t file, int root, int size[3], double h)
{
	int i, err;
	const char *name[] = {
		"vx", "vy", "vz"	};
	for (i = 0; i < 3; i++) {
		err = h5_swrite(file, root, name[i], size[i], h);
		if (err < 0) {
			fprintf(stderr, "%s:%d: h5_swrite failed\n", __FILE__, __LINE__);
			return err;
		}
	}
	return 0;
}

int
h5_close(hid_t file)
{
	return H5Fclose(file);
}

/*
awk '{gsub(/"/, "\\\""); printf "\"%s\\n\"\\\n", $0}' poc/p.xmf
*/
int
h5_xmf(const char *path, int size[3])
{
	FILE *f;
	enum {X, Y, Z};
	int x, y, z, u, v, w;
	char base[MAX_SIZE], full[MAX_SIZE];
	const char s[] = \
"<?xml version=\"1.0\" ?>\n"\
"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"\
"<Xdmf Version=\"2.0\">\n"\
" <Domain>\n"\
"   <Grid GridType=\"Uniform\">\n"\
"     <Time Value=\"1.000000e+02\"/>\n"\
"\n"\
"     <Topology TopologyType=\"3DRectMesh\" Dimensions=\"%d %d %d\"/>\n"\
"\n"\
"     <Geometry GeometryType=\"VxVyVz\">\n"\
"       <DataItem Name=\"mesh_vx\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"\
"        %s.h5:/vx\n"\
"       </DataItem>\n"\
"       <DataItem Name=\"mesh_vy\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"\
"        %s.h5:/vy\n"\
"       </DataItem>\n"\
"       <DataItem Name=\"mesh_vz\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"\
"        %s.h5:/vz\n"\
"       </DataItem>\n"\
"     </Geometry>\n"\
"\n"\
"     <Attribute Name=\"p\" AttributeType=\"Scalar\" Center=\"Cell\">\n"\
"       <DataItem Dimensions=\"%d %d %d 1\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"\
"        p.h5:/data\n"\
"       </DataItem>\n"\
"     </Attribute>\n"\
"   </Grid>\n"\
" </Domain>\n"\
"</Xdmf>\n"\
;
	h5_strcat(path, ".xmf", full);
	h5_basename(path, base);
	f = fopen(full, "w");
	if (f == NULL) {
		fprintf(stderr, "%s:%d: %s : fail to open\n", __FILE__, __LINE__, full);
		return -1;
	}
	x = size[X]; y = size[Y]; z = size[Z];
	u = x + 1; v = y + 1; w = z + 1;
	fprintf(f, s, u, v, w, u, base, v, base, w, base, x, y, z);
	fclose(f);
	return 0;
}
