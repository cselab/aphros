/*
 * Author: Eric R. Schendel <erschend AT ncsu.edu>
 */

#include "compressor_zlib_routines.h"
#include <assert.h>
#include <zlib.h>
#include <stdio.h>

size_t zlib_compress_bound(const size_t _input_buffer_size)
{
    return compressBound(_input_buffer_size);
}

size_t zlib_compress_entire(const size_t _input_buffer_size, char *_input_buffer,
                            const size_t _output_buffer_limit_size, char *_output_buffer,
                            const int _level)
{
    size_t output_buffer_size = 0;

    z_stream strm;

    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;

    assert(deflateInit(&strm, _level) == Z_OK);

    strm.avail_in = _input_buffer_size;
    strm.avail_out = _output_buffer_limit_size;
    strm.next_in = (Bytef*) _input_buffer;
    strm.next_out = (Bytef*) _output_buffer;

    assert(deflate(&strm, Z_FINISH) == Z_STREAM_END);

    output_buffer_size = _output_buffer_limit_size - strm.avail_out;

    assert(output_buffer_size != 0 && output_buffer_size <= _output_buffer_limit_size);

    deflateEnd(&strm);

    return output_buffer_size;
}

void zlib_decompress_entire(const size_t _input_buffer_size, char *_input_buffer,
                            const size_t _output_buffer_size, char *_output_buffer)
{
    z_stream strm;

    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;

    assert(inflateInit(&strm) == Z_OK);

    strm.avail_in = _input_buffer_size;
    strm.avail_out = _output_buffer_size;
    strm.next_in = (Bytef*) _input_buffer;
    strm.next_out = (Bytef*) _output_buffer;

    assert(inflate(&strm, Z_FINISH) == Z_STREAM_END);

    assert(strm.avail_out == 0);

    inflateEnd(&strm);
}

size_t zlib2_compress_entire(const size_t _input_buffer_size, char *_input_buffer,
                            const size_t _output_buffer_size, char *_output_buffer,
                            const int level)
{
    size_t compressed_output_size = _output_buffer_size;
    compress2 ((Bytef*) _output_buffer, &compressed_output_size, (Bytef*) _input_buffer, (uLong) _input_buffer_size, level);
    return compressed_output_size;
}

void zlib2_decompress_entire(const size_t _input_buffer_size, char *_input_buffer,
                            const size_t _output_buffer_size, char *_output_buffer)
{
    uLongf compressed_buffer_size = _output_buffer_size;
    uncompress ((Bytef*) _output_buffer, &compressed_buffer_size, (Bytef*) _input_buffer, compressed_buffer_size);
}
