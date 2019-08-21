int h5_hdf(MPI_Comm, const char *path, int size[3], int start[3], int extent[3], double*);
int h5_xmf(const char *path, const char *name, double origin[3], double spacing, int size[3]);
int h5_silence(void);
