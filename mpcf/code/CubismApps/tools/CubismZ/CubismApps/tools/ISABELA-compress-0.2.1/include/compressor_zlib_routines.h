/*
 * Author: Eric R. Schendel <erschend AT ncsu.edu>
 */

#ifndef COMPRESSOR_ZLIB_ROUTINES_H
#define	COMPRESSOR_ZLIB_ROUTINES_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <stddef.h>

/***
 * Returns the byte size necessary for a zlib output buffer to compress the
 * worst case scenario for a given input buffer size.
 */
size_t zlib_compress_bound(const size_t input_buffer_size);

/***
 * Returns output buffer size.
 *
 * Compression operation is done on the entire buffer as single operation
 * with the assumption that the output buffer memory allocation will be large
 * enough (see, zlib_compress_bound).
 */
size_t zlib_compress_entire(const size_t input_buffer_size, char *input_buffer,
                            const size_t output_buffer_limit_size, char *output_buffer,
                            const int level);

/***
 * Decompression operation is done on the entire buffer as single operation
 * where the caller must know the input and output buffer size before hand.
 * If the decompressed output buffer size has not been previously determined,
 * then this function is useless.
 */
void zlib_decompress_entire(const size_t input_buffer_size, char *input_buffer,
                            const size_t output_buffer_size, char *output_buffer);

size_t zlib2_compress_entire(const size_t _input_buffer_size, char *_input_buffer,
                            const size_t _output_buffer_size, char *_output_buffer, 
                            const int level);

void zlib2_decompress_entire(const size_t _input_buffer_size, char *_input_buffer,
                            const size_t _output_buffer_size, char *_output_buffer);
#ifdef	__cplusplus
}
#endif

#endif	/* COMPRESSOR_ZLIB_ROUTINES_H */
