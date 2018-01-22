/* minimal code example showing how to call the zfp (de)compressor */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "myisa.h"

//#define MAX(a,b) (a)>(b)?(a):(b)
double isa_rate = 0.0;

int main(int argc, char* argv[])
{
  if (argc == 2) isa_rate = atof(argv[1]);
  /* allocate 100x100x100 array of doubles or so */
  int nx = 32;
  int ny = 32;
  int nz = 32;

  int is_float = 1;
  size_t init_bytes = nx*ny*nz*sizeof(float);
  float* array = (float *)malloc(init_bytes);
  float* array2 = (float *)malloc(init_bytes);

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

  /* compress or decompress array */
  unsigned char * buffer = (unsigned char *)malloc(nx * ny * nz * (is_float? sizeof(float): sizeof(double)));
  int zbytes = 0;

//static int isa_compress_buffer(char *in_buffer, int nx, int ny, int nz, int rate, int is_float, char *out_buffer, int *zbytes)

  int status = isa_compress_buffer((char *)array, nx, ny, nz, isa_rate, is_float, (char *)buffer, &zbytes);
  printf("status = %d, bytes: from %ld to %ld  (CR=%.2f)\n", status, init_bytes, zbytes, 1.0*init_bytes/zbytes);

  int ubytes = 0;
  int status2 = isa_decompress_buffer((char *)array2, nx, ny, nz, isa_rate, is_float, (char *)buffer, zbytes, &ubytes);
  printf("status2 = %d, bytes: from %ld to %ld\n", status2, zbytes, ubytes);

  double max_err = 0;
  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++) {
	int idx = i + nx * (j + ny * k);
	double diff = fabs(array[idx]-array2[idx]);
	double maxv = MAX(fabs(array[idx]), fabs(array2[idx]));
	double err = 100.0*diff/maxv;
	/*printf("%d: %f vs %f : %.2f%%\n", idx, array[idx], array2[idx], err);*/
	if (err > max_err) max_err = err;
      }
  printf("ACC=%f CR=%.2f max_err = %.2f%%\n", isa_rate, 1.0*init_bytes/zbytes, max_err);

  return 0;

}
