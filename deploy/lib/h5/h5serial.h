int h5_read_xmf(const char *path, /**/ char *name, double origin[3], double *spacing, int size[3]);
int h5_read_hdf(const char *path, /**/ int size[3], double**);
int h5_serial_hdf(const char *path, const int size[3], const double*);
