/*
** fpzip version 1.1.0, June 8, 2014
** Part of the LOCAL Toolkit, UCRL-CODE-232243
** Written by Peter Lindstrom, Lawrence Livermore National Laboratory
**
**
**                               SUMMARY
**
** fpzip is a library for lossless or lossy compression of 2D and 3D arrays
** of single- or double-precision floating-point scalars.  Although linear
** (1D) streams of scalars may be compressed as well, the library has
** primarily been designed to exploit higher-dimensional structure in the
** data, and may not perform well on 1D data or on data not stored as a
** spatially coherent contiguous array, with 'smoothly' varying values.
**
** The library supports compression to either a file or to main memory, and
** allows specifying how many bits of precision to retain by truncating
** each floating-point value and discarding the least significant bits; the
** remaining bits are compressed losslessly.  The precision is limited to
** integers 2-32 for floats.  For doubles, precisions 4-64 are supported in
** increments of two bits.  The decompressed data is returned in full
** precision, with any truncated bits zeroed.
**
** Because floating-point arithmetic may be affected by factors such as
** register precision, rounding mode, and compiler optimizations,
** precautions have been taken to ensure correctness and portability via a
** set of compile-time macros.  For example, it is possible to specify that
** floating-point operations be emulated via integer arithmetic, or to
** treat the binary representation of floating-point numbers as integers.
** Please consult the Makefile for choosing among these settings.  The
** compressor works correctly on the IEEE 754 floating-point format, though
** no particular assumption is made on the floating-point representation
** other than the most significant bit being the sign bit.  Special values
** such as infinities, NaNs, and denormalized numbers should be handled
** correctly by the compressor in lossless mode.
**
** For convenience, functions are provided for reading and writing a header
** that stores version specific information and meta data such as array
** dimensions.  Such meta data must be provided to both compressor and
** decompressor.  It is up to the caller to decide whether or not to use
** such fpzip headers, as this information may already be available
** externally from container formats like HDF5.  If fpzip headers are not
** used, then the FPZ fields must be set by the caller in both read and
** write mode.
**
** A single compressed stream may store multiple contiguous fields (e.g.
** for multiple arrays with the same dimensions that represent different
** variables).  Similarly, a stream may store multiple arrays of different
** dimensions and types, possibly with one header per array.  It is up to
** the caller to interleave read/write calls that perform (de)compression
** of floating-point data with read/write calls of header data.
**
** The return value of each function should be checked in case invalid
** arguments are passed or a run-time error occurs.  In this case, the
** variable fpzip_errno is set and can be examined to determine the cause
** of the error.
**
** fpzip was developed as part of the LOCAL LDRD project at LLNL, and may
** be freely used and distributed for noncommercial purposes.  The core
** library is written in C++ and applications need to be linked with a C++
** linker.  The library can, however, be called from C.  For further
** information and bug reports, please e-mail pl@llnl.gov.
**
**
**                                NOTICE
**
** This work was produced at the Lawrence Livermore National Laboratory
** (LLNL) under contract no. DE-AC-52-07NA27344 (Contract 44) between
** the U.S. Department of Energy (DOE) and Lawrence Livermore National
** Security, LLC (LLNS) for the operation of LLNL.  The rights of the
** Federal Government are reserved under Contract 44.
**
**
**                              DISCLAIMER
**
** This work was prepared as an account of work sponsored by an agency of
** the United States government.  Neither the United States government nor
** Lawrence Livermore National Security, LLC, nor any of their employees
** makes any warranty, expressed or implied, or assumes any legal liability
** or responsibility for the accuracy, completeness, or usefulness of any
** information, apparatus, product, or process disclosed, or represents
** that its use would not infringe privately owned rights.  Reference
** herein to any specific commercial product, process, or service by trade
** name, trademark, manufacturer, or otherwise does not necessarily
** constitute or imply its endorsement, recommendation, or favoring by the
** United States government or Lawrence Livermore National Security, LLC.
** The views and opinions of authors expressed herein do not necessarily
** state or reflect those of the United States government or Lawrence
** Livermore National Security, LLC, and shall not be used for advertising
** or product endorsement purposes.
**
**
**                            COMMERCIAL USE
**
** Commercialization of this product is prohibited without notifying the
** Department of Energy (DOE) or Lawrence Livermore National Security.
*/

#ifndef FPZIP_H
#define FPZIP_H

#include <stdio.h>

#define FPZIP_TYPE_FLOAT  0 // single-precision data (see FPZ.type)
#define FPZIP_TYPE_DOUBLE 1 // double-precision data

#ifdef __cplusplus
extern "C" {
#endif

/* array meta data and stream handle */
typedef struct {
  int type; /* single (0) or double (1) precision */
  int prec; /* number of bits of precision (zero = full) */
  int nx;   /* number of x samples */
  int ny;   /* number of y samples */
  int nz;   /* number of z samples */
  int nf;   /* number of fields */
} FPZ;

/* associate file with compressed input stream */
FPZ*                  /* compressed stream */
fpzip_read_from_file(
  FILE* file          /* binary input stream */
);

/* associate memory buffer with compressed input stream */
FPZ*                  /* compressed stream */
fpzip_read_from_buffer(
  const void* buffer  /* pointer to compressed input data */
);

/* read FPZ meta data (use only if previously written) */
int                   /* nonzero upon success */
fpzip_read_header(
  FPZ* fpz            /* compressed stream */
);

/* decompress array */
size_t                /* number of compressed bytes read (zero = error) */
fpzip_read(
  FPZ*  fpz,          /* compressed stream */
  void* data          /* uncompressed floating-point data */
);

/* close input stream and deallocate fpz */
void
fpzip_read_close(
  FPZ* fpz            /* compressed stream */
);

/* associate file with compressed output stream */
FPZ*                  /* compressed stream */
fpzip_write_to_file(
  FILE* file          /* binary output stream */
);

/* associate memory buffer with compressed output stream */
FPZ*                  /* compressed stream */
fpzip_write_to_buffer(
  void*  buffer,      /* pointer to compressed output data */
  size_t size         /* size of allocated storage for buffer */
);

/* write FPZ meta data */
int                   /* nonzero upon success */
fpzip_write_header(
  FPZ* fpz            /* compressed stream */
);

/* compress array */
size_t                /* number of compressed bytes written (zero = error) */
fpzip_write(
  FPZ*        fpz,    /* compressed stream */
  const void* data    /* uncompressed floating-point data */
);

/* close output stream and deallocate fpz */
void
fpzip_write_close(
  FPZ* fpz            /* compressed stream */
);

/*
** Error codes.
*/

typedef enum {
  fpzipSuccess             = 0, /* no error */
  fpzipErrorReadStream     = 1, /* cannot read stream */
  fpzipErrorWriteStream    = 2, /* cannot write stream */
  fpzipErrorBadFormat      = 3, /* magic mismatch; not an fpz stream */
  fpzipErrorBadVersion     = 4, /* fpz format version not supported */
  fpzipErrorBadPrecision   = 5, /* precision not supported */
  fpzipErrorBufferOverflow = 6  /* compressed buffer overflow */
} fpzipError;

extern fpzipError fpzip_errno;     /* error code */
extern const char* fpzip_errstr[]; /* error message indexed by fpzip_errno */

#ifdef __cplusplus
}
#endif

#endif
