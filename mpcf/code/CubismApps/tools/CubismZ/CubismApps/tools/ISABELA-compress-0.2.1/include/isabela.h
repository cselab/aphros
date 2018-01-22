/*
 * Author: Sriram Lakshminarasimhan <slakshm2@ncsu.edu>
 */

#ifndef __ISABELA_HEADER__
#define __ISABELA_HEADER__

#include "utils.h"
#include <stdint.h>

/*
 * The following ISABELA API is based on the zlib library API.
 */

enum ISABELA_status
{
    ISABELA_SUCCESS = 0, /* when a function returns ISABELA_SUCCESS, everything is OK */
    ISABELA_ERROR   = 1  /* this is a catch all for ISABELA function errors */
};

enum ISABELA_flush
{
    ISABELA_FINISH = 0   /* only flush option available so far (similar to zlib library) */
};

enum ISABELA_type
{
    ISABELA_BSPLINES = 0,   /* Approximate sorted windows using BSplines */
    ISABELA_WAVELETS = 1    /* Approximate sorted windows using Haar Wavelets */
};

struct isabela_config
{
    uint32_t window_size;       /*(= W_0 in paper). Number of floating-point elements in each window. 
                                  Default window size is 1024 elements. */
    uint32_t ncoefficients;     /*(= C in paper). Number of co-efficients to use to approximate the 
                                  sorted values within each window. Default number of coefficients is 
                                  30. */
    float error_rate;           /*(= epsilon in paper). Max relative error between the approximated 
                                  and original of each point. Use error_rate = 0 to ignore error 
                                  encoding. Default error rate is 1%. */
    uchar element_byte_size;    /* 4 bytes for single-precision float, 8 bytes for double-precision */
    uchar transform;            /* Switch between wavelets and BSplines. */
};

struct isabela_stream
{
    void *next_in;          /* next input byte */
    uint32_t avail_in;  /* number of bytes available at next_in */

    void *next_out;         /* next output byte */
    uint32_t avail_out; /* number of bytes available at next_out */

    /* The following pointers reference the sections within next_out buffer */
    /* order of buffer is meta_ptr, incomp_ptr, comp_ptr, then end_ptr      */
    /* note: this interim support until ISABELA_flush support provided     */

    void *meta_ptr;         /* pointer to meta data in buffer (first byte/segment of next_out) */
    void *index_ptr;        /* pointer to index data (next byte after meta segment) */
    void *coefficients_ptr; /* pointer to coefficients data (next byte after index segment) */
    void *error_ptr;        /* pointer to compressed error data (next byte after coefficients segment) */
    void *error_size_ptr;   /* pointer to compressed error data (next byte after coefficients segment) */
    void *end_ptr;          /* pointer to the end of the data (next byte after error segment) */

    struct isabela_config *config;  /* pointer to ISABELA configuration object*/
};

/*
 * isabelaDeflateInit is required for initializing isabela_stream data structure before any other
 * ISABELA deflate functions are called.
 *
 * Input:
 *   ib_stream:     pointer to isabela_stream structure required by all ISABELA Deflate functions.
 *   datatype_size: byte size of datatype representing each element in a dataset (ie, double=8).
 *   ib_preference: ISABELA compression preference (ie, ISABELA_SPEED).
 */
enum ISABELA_status
isabelaDeflateInit (struct isabela_stream *ib_stream,
                    uint32_t datatype_size,
                    struct isabela_config *is_config);

/*
 * isabelaDeflate applies ISABELA lossy compression dataset preconditioner and compresses
 * input buffer (see 'avail_in' above) and generates output to 'avail_out'.
 *
 * Input:
 *   ib_stream: pointer to isabela_stream structure initialized by isabelaDeflateInit.
 *   ib_flush:  currently only supports ISABELA_FINISH (similar to zlib Z_FINISH).
 */
enum ISABELA_status
isabelaDeflate (struct isabela_stream *ib_stream,
                enum ISABELA_flush ib_flush);

/*
 * isabelaDeflateEnd deallocates isabela_stream data structure when no longer needed.
 *
 * Input:
 *   i_stream: pointer to isabela_stream structure initialized by isabelaDeflateInit.
 */
enum ISABELA_status
isabelaDeflateEnd (struct isabela_stream *ib_stream);

/*
 * isabelaInflateInit is required for initializing isabela_stream data structure before any other
 * ISABELA inflate functions are called.
 *
 * Input:
 *   ib_stream:     pointer to isabela_stream structure required by all ISABELA Inflate functions.
 *   datatype_size: byte size of datatype representing each element in a dataset (ie, double=8).
 */
enum ISABELA_status
isabelaInflateInit (struct isabela_stream *ib_stream,
                    uint32_t datatype_size,
                    struct isabela_config *is_config);

/*
 * isabelaInflate reverses lossy compression applied by isabelaDeflate using the
 * input buffer (see 'avail_in' above) and generates output to 'avail_out'.
 *
 * Input:
 *   ib_stream: pointer to isabela_stream structure initialized by isabelaInflateInit.
 *   ib_flush:  currently only supports ISABELA_FINISH (similar to zlib Z_FINISH).
 */
enum ISABELA_status
isabelaInflate (struct isabela_stream *ib_stream,
                enum ISABELA_flush ib_flush);

/*
 * isabelaInflateEnd deallocates isabela_stream data structure when no longer needed.
 *
 * Input:
 *   i_stream: pointer to isabela_stream structure initialized by isabelaInflateInit.
 */
enum ISABELA_status
isabelaInflateEnd (struct isabela_stream *ib_stream);

enum ISABELA_status isabelaPrecisionTest (double *original, double *approximated, uint32_t n, float error_rate);

#endif /* __ISABELA_HEADER__ */
