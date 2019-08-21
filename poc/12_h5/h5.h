int h5_close(hid_t);
hid_t h5_open(MPI_Comm, const char*);
int h5_data(hid_t, int size[3], int start[3], int extent[3], double*);
int h5_xmf(const char *path, const char *name, double origin[3], double spacing, int[3]);
