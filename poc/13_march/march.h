#ifdef __cplusplus
extern "C" {
#endif

int march_cube(double u[8], /**/ int *ntri, double tri[4 * 3 * 3]);
int march_tetrahedron(double u[8], /**/ int *ntri, double tri[4 * 3 * 3]);

#ifdef __cplusplus
}
#endif
