int h5_hdf(MPI_Comm, const char *path, const int size[3], const int start[3], const int extent[3], const double*);
int h5_xmf(const char *path, const char *name, const double origin[3], double spacing, const int size[3]);
int h5_silence(void);
