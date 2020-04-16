int h5_xmf(
    const char* path, const char* name, const double origin[3], double spacing,
    const int size[3]);
int h5_silence(void);
int h5_serial_hdf(const char* path, const int size[3], const double* data);
int h5_read_hdf(const char* path, /**/ int size[3], double**);
int h5_read_xmf(
    const char*, /**/ char* name, double origin[3], double spacing[3],
    int size[3]);
