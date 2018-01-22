/*
 * Author: Sriram Lakshminarasimhan <slakshm2 AT ncsu.edu>
 */

#include "isabela.h"

#include "compressor_zlib_routines.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "error.h"
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "bitstream.h"
#include <assert.h>

#include "spline_fit.h"
#include "wavelets.h"

void print_isabela_config (struct isabela_config *is_config)
{
    static int once = 0;
    if (once == 1) return;
    else once = 1; 
    printf ("\n");
    printf ("Window Size : %d\n", is_config->window_size);
    printf ("Coefficients: %d\n", is_config->ncoefficients);
    printf ("Error Rate  : %f\n", is_config->error_rate);
    printf ("\n");
}

/* PRIVATE FUNCTIONS */

enum ISABELA_status
isabela_init(struct isabela_stream* is_stream,
            uint32_t datatype_size,
            struct isabela_config* is_config)
{

    if (is_stream == NULL) {
        return ISABELA_ERROR;
    }

    if (datatype_size == 0 || datatype_size > 255) {
        return ISABELA_ERROR;
    }

    is_stream->next_in = NULL;
    is_stream->avail_in = 0;
    is_stream->next_out = NULL;
    is_stream->avail_out = 0;

    is_stream->config = (struct isabela_config *) malloc (sizeof (struct isabela_config));

    if (is_stream->config == NULL) {
        return ISABELA_ERROR;
    }

    // Default values when config object is not passed
    if (is_config == NULL) {
        is_stream->config->window_size = 1024;
        is_stream->config->ncoefficients = 30;
        is_stream->config->error_rate = 1;
        is_stream->config->element_byte_size = sizeof (double);
        is_stream->config->transform = ISABELA_BSPLINES;
    } else {
        is_stream->config->window_size = is_config->window_size;
        is_stream->config->ncoefficients = is_config->ncoefficients;
        is_stream->config->error_rate = is_config->error_rate;
        is_stream->config->element_byte_size = is_config->element_byte_size;
        is_stream->config->transform = is_config->transform;
    }

    return ISABELA_SUCCESS;
}


enum ISABELA_status
isabela_end(struct isabela_stream* is_stream)
{
    free(is_stream->config);

    return ISABELA_SUCCESS;
}

uint32_t isabela_get_config_size (struct isabela_config *config)
{
    return sizeof (config->window_size)             +
            sizeof (config->ncoefficients)          +
            sizeof (config->element_byte_size)      +
            sizeof (config->error_rate)             +
            sizeof (config->transform);
}

void isabela_serialize_config (struct isabela_config *config, uchar *out_buffer_uc)
{
    // Set and advance
    SET (uint32_t, out_buffer_uc, config->window_size);
    out_buffer_uc += sizeof (uint32_t);

    SET (uint32_t, out_buffer_uc, config->ncoefficients);
    out_buffer_uc += sizeof (uint32_t);

    SET (uchar, out_buffer_uc, config->element_byte_size);
    out_buffer_uc += sizeof (uchar);

    SET (float, out_buffer_uc, config->error_rate);
    out_buffer_uc += sizeof (float);

    SET (uchar, out_buffer_uc, config->transform);
    out_buffer_uc += sizeof (float);
}

void isabela_deserialize_config (struct isabela_config *config, uchar *out_buffer_uc)
{
    // Set and advance
    GET (uint32_t, out_buffer_uc, config->window_size);
    out_buffer_uc += sizeof (uint32_t);

    GET (uint32_t, out_buffer_uc, config->ncoefficients);
    out_buffer_uc += sizeof (uint32_t);

    GET (uchar, out_buffer_uc, config->element_byte_size);
    out_buffer_uc += sizeof (uchar);

    GET (float, out_buffer_uc, config->error_rate);
    out_buffer_uc += sizeof (float);

    GET (uchar, out_buffer_uc, config->transform);
    out_buffer_uc += sizeof (float);
}

#include "isabela_float.c"
#include "isabela_double.c"

/* PUBLIC FUNCTIONS */

enum ISABELA_status
isabelaDeflateInit(struct isabela_stream* is_stream,
                  uint32_t datatype_size,
                  struct isabela_config *is_config)
{
    enum ISABELA_status status = isabela_init(is_stream, datatype_size, is_config);

    return status;
}

enum ISABELA_status
isabelaDeflateEnd(struct isabela_stream* ib_stream)
{
    return isabela_end(ib_stream);
}

enum ISABELA_status
isabelaInflateInit(struct isabela_stream* ib_stream,
                  uint32_t datatype_size,
                  struct isabela_config *is_config)
{
    return isabela_init(ib_stream, datatype_size, is_config);
}

enum ISABELA_status
isabelaInflateEnd(struct isabela_stream* ib_stream)
{
    return isabela_end(ib_stream);
}

enum ISABELA_status
isabelaInflate(struct isabela_stream* ib_stream, enum ISABELA_flush ib_flush)
{
    if (ib_stream->config->element_byte_size == sizeof (double)) {
        return isabelaInflateDouble (ib_stream, ib_flush);
    } else if (ib_stream->config->element_byte_size == sizeof (float)) {
        return isabelaInflateFloat (ib_stream, ib_flush);
    }
}

enum ISABELA_status
isabelaDeflate(struct isabela_stream* ib_stream, enum ISABELA_flush ib_flush)
{
    if (ib_stream->config->element_byte_size == sizeof (double)) {
        return isabelaDeflateDouble (ib_stream, ib_flush);
    } else if (ib_stream->config->element_byte_size == sizeof (float)) {
        return isabelaDeflateFloat (ib_stream, ib_flush);
    }
}
