#ifdef __cplusplus
extern "C" {
#endif

enum { MARCH_NTRI = 48 };
int march_cube(double u[8], /**/ int *ntri, double *tri);
int march_cube_location(/**/ int*, int*, double*);

#ifdef __cplusplus
}
#endif
