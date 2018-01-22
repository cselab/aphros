/* minimal code example showing how to call the zfp (de)compressor */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "myzfp.h"

double ZFP_ACC = 0.0;

int main(int argc, char* argv[])
{
  if (argc == 2) ZFP_ACC = atof(argv[1]);
  /* allocate 100x100x100 array of doubles or so */
  int nx = 32;
  int ny = 32;
  int nz = 32;

#if 0
  int is_float = 0;
  size_t init_bytes = nx*ny*nz*sizeof(double);
  double* array = malloc(init_bytes);
  double* array2 = malloc(init_bytes);

  /* initialize array to be compressed */
  int i, j, k;
  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++) {
        double x = 2.0 * i / nx;
        double y = 2.0 * j / ny;
        double z = 2.0 * k / nz;
        array[i + nx * (j + ny * k)] = exp(-(x * x + y * y + z * z));
      }
#else
  int is_float = 1;
  size_t init_bytes = nx*ny*nz*sizeof(float);
  float* array = malloc(init_bytes);
  float* array2 = malloc(init_bytes);

  /* initialize array to be compressed */
  int i, j, k;
  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++) {
        float x = 2.0 * i / nx;
        float y = 2.0 * j / ny;
        float z = 2.0 * k / nz;
        array[i + nx * (j + ny * k)] = exp(-(x * x + y * y + z * z));
      }
#endif

  /* compress or decompress array */
  unsigned char * buffer = malloc(nx * ny * nz * (is_float? sizeof(float): sizeof(double)));
  size_t zbytes = 0;

  int status = zfp_compress_buffer(array, nx, ny, nz, ZFP_ACC, is_float, buffer, &zbytes);
  printf("status = %d, bytes: from %ld to %ld  (CR=%.2f)\n", status, init_bytes, zbytes, 1.0*init_bytes/zbytes);

  size_t ubytes = 0;
  int status2 = zfp_decompress_buffer(array2, nx, ny, nz, ZFP_ACC, is_float, buffer, zbytes, &ubytes);
  printf("status2 = %d, bytes: from %ld to %ld\n", status2, zbytes, ubytes);

  double max_err = 0;
  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++) {
	int idx = i + nx * (j + ny * k);
	double diff = fabs(array[idx]-array2[idx]);
	double maxv = max(fabs(array[idx]), fabs(array2[idx]));
	double err = 100.0*diff/maxv;
	/*printf("%d: %f vs %f : %.2f%%\n", idx, array[idx], array2[idx], err);*/
	if (err > max_err) max_err = err;
      }
  printf("ACC=%f CR=%.2f max_err = %.2f%%\n", ZFP_ACC, 1.0*init_bytes/zbytes, max_err);

  return 0;

}
