#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "fpzip.h"

// return current time
static double
now()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1e-6 * tv.tv_usec;
}

// compress 'inbytes' bytes from 'data' to stream 'FPZ'
static int
compress(FPZ* fpz, const void* data, size_t inbytes, const char* medium)
{
  fprintf(stderr, "compressing to %s\n", medium);
  double t = now();
#ifndef WITHOUT_HEADER
  // write header (optional)
  if (!fpzip_write_header(fpz)) {
    fprintf(stderr, "cannot write header: %s\n", fpzip_errstr[fpzip_errno]);
    return 0;
  }
#endif
  // perform actual compression
  size_t outbytes = fpzip_write(fpz, data);
  if (!outbytes) {
    fprintf(stderr, "compression failed: %s\n", fpzip_errstr[fpzip_errno]);
    return 0;
  }
  t = now() - t;
  fprintf(stderr, "in=%zu out=%zu ratio=%.2f seconds=%.3f MB/s=%.3f\n", inbytes, outbytes, (double)inbytes / outbytes, t, (double)inbytes / (1024 * 1024 * t));
  return 1;
}

// decompress 'inbytes' bytes to 'data' from stream 'FPZ'
static int
decompress(FPZ* fpz, void* data, size_t inbytes, const char* medium)
{
  fprintf(stderr, "decompressing from %s\n", medium);
  double t = now();
#ifndef WITHOUT_HEADER
  // read header (if previously written)
  if (!fpzip_read_header(fpz)) {
    fprintf(stderr, "cannot read header: %s\n", fpzip_errstr[fpzip_errno]);
    return 0;
  }
  // make sure array size stored in header matches expectations
  size_t size = (fpz->type == FPZIP_TYPE_FLOAT ? sizeof(float) : sizeof(double));
  if (size * fpz->nx * fpz->ny * fpz->nz * fpz->nf != inbytes) {
    fprintf(stderr, "array size does not match dimensions from header\n");
    return 0;
  }
#endif
  // perform actual decompression
  if (!fpzip_read(fpz, data)) {
    fprintf(stderr, "decompression failed: %s\n", fpzip_errstr[fpzip_errno]);
    return EXIT_FAILURE;
  }
  t = now() - t;
  fprintf(stderr, "seconds=%.3f MB/s=%.3f\n", t, (double)inbytes / (1024 * 1024 * t));
  return 1;
}

// validate float data
static int
validate_float(const float* data, const float* copy, size_t count)
{
  fprintf(stderr, "validating\n");
  for (size_t i = 0; i < count; i++)
    if (data[i] != copy[i]) {
      fprintf(stderr, "validation failed at i=%zu\n", i);
      return 0;
    }
  return 1;
}

// validate double data
static int
validate_double(const double* data, const double* copy, size_t count)
{
  fprintf(stderr, "validating\n");
  for (size_t i = 0; i < count; i++)
    if (data[i] != copy[i]) {
      fprintf(stderr, "validation failed at i=%zu\n", i);
      return 0;
    }
  return 1;
}

// compute error for float data
static double
error_float(const float* data, const float* copy, size_t count)
{
  fprintf(stderr, "computing error\n");
  double e = 0;
  for (size_t i = 0; i < count; i++)
    e += (data[i] - copy[i]) * (data[i] - copy[i]);
  return sqrt(e / count);
}

// compute error for double data
static double
error_double(const double* data, const double* copy, size_t count)
{
  fprintf(stderr, "computing error\n");
  double e = 0;
  for (size_t i = 0; i < count; i++)
    e += (data[i] - copy[i]) * (data[i] - copy[i]);
  return sqrt(e / count);
}

// verify lossless compression or measure loss
static int
validate(const void* data, const void* copy, size_t count, int type, int lossless)
{
  if (lossless) {
    if (!(type == FPZIP_TYPE_FLOAT ? validate_float((const float*)data, (const float*)copy, count) : validate_double((const double*)data, (const double*)copy, count)))
      return 0;
  }
  else
    fprintf(stderr, "rmse=%g\n", (type == FPZIP_TYPE_FLOAT ? error_float((const float*)data, (const float*)copy, count) : error_double((const double*)data, (const double*)copy, count)));
  return 1;
}

int main(int argc, char* argv[])
{
  int type = FPZIP_TYPE_FLOAT;
  int prec = 0;
  int nx = 1;
  int ny = 1;
  int nz = 1;
  int nf = 1;
  char* infile = 0;
  char* outfile = 0;
  char* reconfile = 0;

  // parse command-line options
  switch (argc) {
    case 10:
      reconfile = argv[9];
      /*FALLTHROUGH*/
    case 9:
      outfile = argv[8];
      /*FALLTHROUGH*/
    case 8:
      infile = argv[7];
      /*FALLTHROUGH*/
    case 7:
      if (sscanf(argv[6], "%d", &prec) != 1)
        goto usage;
      /*FALLTHROUGH*/
    case 6:
      if (sscanf(argv[5], "%d", &nf) != 1)
        goto usage;
      /*FALLTHROUGH*/
    case 5:
      if (sscanf(argv[4], "%d", &nz) != 1)
        goto usage;
      /*FALLTHROUGH*/
    case 4:
      if (sscanf(argv[3], "%d", &ny) != 1)
        goto usage;
      /*FALLTHROUGH*/
    case 3:
      if (sscanf(argv[2], "%d", &nx) != 1)
        goto usage;
      /*FALLTHROUGH*/
    case 2:
      if (!strcmp(argv[1], "-f"))
        type = FPZIP_TYPE_FLOAT;
      else if (!strcmp(argv[1], "-d"))
        type = FPZIP_TYPE_DOUBLE;
      else
        goto usage;
      break;
    default:
    usage:
      fprintf(stderr, "Usage: testfpzip <-f|-d> [nx [ny [nz [nf [prec [infile [outfile [reconfile]]]]]]]]\n");
      return EXIT_FAILURE;
  }

  // initialize
  fprintf(stderr, "testing '%s a[%d][%d][%d][%d]'\n", type == FPZIP_TYPE_FLOAT ? "float" : "double", nf, nz, ny, nx);
  size_t count = (size_t)nx * ny * nz * nf;
  size_t size = (type == FPZIP_TYPE_FLOAT ? sizeof(float) : sizeof(double));
  size_t inbytes = count * size;
  size_t bufbytes = 1024 + inbytes;
  if (prec == 0)
    prec = CHAR_BIT * size;
  int lossless = (prec == CHAR_BIT * size);

  // allocate buffers
#ifdef __cplusplus
  void* data = (type == FPZIP_TYPE_FLOAT ? static_cast<void*>(new float[count]) : static_cast<void*>(new double[count]));
  void* copy = (type == FPZIP_TYPE_FLOAT ? static_cast<void*>(new float[count]) : static_cast<void*>(new double[count]));
  void* buffer = static_cast<void*>(new unsigned char[bufbytes]);
#else
  void* data = malloc(inbytes);
  void* copy = malloc(inbytes);
  void* buffer = malloc(bufbytes);
#endif

  // initialize scalar field
  if (infile) {
    // read raw data
    fprintf(stderr, "reading input file\n");
    FILE* file = fopen(infile, "rb");
    if (!file) {
      fprintf(stderr, "cannot open input file\n");
      return EXIT_FAILURE;
    }
    if (fread(data, size, count, file) != count) {
      fprintf(stderr, "read failed\n");
      return EXIT_FAILURE;
    }
    fclose(file);
  }
  else {
    // set scalar field to product of cosines
    const double pi = acos(-1.0);
    if (type == FPZIP_TYPE_FLOAT) {
      float* p = (float*)data;
      for (int f = 0; f < nf; f++)
        for (int z = 0; z < nz; z++)
          for (int y = 0; y < ny; y++)
            for (int x = 0; x < nx; x++)
              *p++ = cos(2 * pi * x / nx) * cos(2 * pi * y / ny) * cos(2 * pi * z / nz);
    }
    else {
      double* p = (double*)data;
      for (int f = 0; f < nf; f++)
        for (int z = 0; z < nz; z++)
          for (int y = 0; y < ny; y++)
            for (int x = 0; x < nx; x++)
              *p++ = cos(2 * pi * x / nx) * cos(2 * pi * y / ny) * cos(2 * pi * z / nz);
    }
  }

  // compress to memory
  FPZ* fpz = fpzip_write_to_buffer(buffer, bufbytes);
  fpz->type = type;
  fpz->prec = prec;
  fpz->nx = nx;
  fpz->ny = ny;
  fpz->nz = nz;
  fpz->nf = nf;
  if (!compress(fpz, data, inbytes, "memory"))
    return EXIT_FAILURE;
  fpzip_write_close(fpz);

  // decompress from memory
  fpz = fpzip_read_from_buffer(buffer);
#ifdef WITHOUT_HEADER
  // manually set array parameters since header is not stored
  fpz->type = type;
  fpz->prec = prec;
  fpz->nx = nx;
  fpz->ny = ny;
  fpz->nz = nz;
  fpz->nf = nf;
#endif
  if (!decompress(fpz, copy, inbytes, "memory"))
    return EXIT_FAILURE;
  fpzip_read_close(fpz);

  // perform validation
  if (!validate(data, copy, count, type, lossless))
    return EXIT_FAILURE;

  if (outfile) {
    // compress to file
    FILE* file = fopen(outfile, "wb");
    if (!file) {
      fprintf(stderr, "cannot create compressed file\n");
      return EXIT_FAILURE;
    }
    fpz = fpzip_write_to_file(file);
    fpz->type = type;
    fpz->prec = prec;
    fpz->nx = nx;
    fpz->ny = ny;
    fpz->nz = nz;
    fpz->nf = nf;
    if (!compress(fpz, data, inbytes, "file"))
      return EXIT_FAILURE;
    fpzip_write_close(fpz);
    fclose(file);

    // decompress from file
    file = fopen(outfile, "rb");
    if (!file) {
      fprintf(stderr, "cannot open compressed file\n");
      return EXIT_FAILURE;
    }
    fpz = fpzip_read_from_file(file);
#ifdef WITHOUT_HEADER
    // manually set array parameters since header is not stored
    fpz->type = type;
    fpz->prec = prec;
    fpz->nx = nx;
    fpz->ny = ny;
    fpz->nz = nz;
    fpz->nf = nf;
#endif
    if (!decompress(fpz, copy, inbytes, "file"))
      return EXIT_FAILURE;
    fpzip_read_close(fpz);
    fclose(file);

    // perform validation
    if (!validate(data, copy, count, type, lossless))
      return EXIT_FAILURE;

    if (reconfile) {
      // write reconstructed data to file
      FILE* file = fopen(reconfile, "wb");
      if (!file) {
        fprintf(stderr, "cannot open reconstructed file\n");
        return EXIT_FAILURE;
      }
      if (fwrite(copy, size, count, file) != count) {
        fprintf(stderr, "write failed\n");
        return EXIT_FAILURE;
      }
      fclose(file);
    }
  }

  // deallocate buffers
#ifdef __cplusplus
  if (type == FPZIP_TYPE_FLOAT) {
    delete[] static_cast<float*>(data);
    delete[] static_cast<float*>(copy);
  }
  else {
    delete[] static_cast<double*>(data);
    delete[] static_cast<double*>(copy);
  }
  delete[] static_cast<unsigned char*>(buffer);
#else
  free(data);
  free(copy);
  free(buffer);
#endif

  fprintf(stderr, "OK\n");

  return 0;
}
