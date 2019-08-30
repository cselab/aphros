#ifdef __cplusplus
extern "C" {
#endif
    enum { MARCH_NTRI = 48 };
    int march_cube(double u[8], /**/ int *ntri, double *tri);
    int march_tetrahedron(double u[8], /**/ int *ntri, double *tri);
#ifdef __cplusplus
}
#endif
