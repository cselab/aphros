// Created by Petr Karnakov on 31.01.2021
// Copyright 2021 ETH Zurich

#ifdef __cplusplus
extern "C" {
#endif

enum { MARCH_NTRI = 48 };
extern const int MARCH_O[][3];
int march_cube(double u[8], /**/ int* ntri, double* tri);
int march_cube_location(
    double u[8], /**/ int* ntri, double* tri, int* x, int* y, double* a);
#ifdef __cplusplus
}
#endif
