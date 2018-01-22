/*
 * Author: John Jenkins <jpjenki2 AT ncsu.edu>
 */
#ifndef BITSTREAM_H
#define BITSTREAM_H

#include <stdint.h>

//fetch nitems number of elements, each with bits_per_item bits, from packed_buffer to
// unpacked buffer, starting at offset in packed buffer. offset is set to the next uint32_t
// to be read from the packed buffer
void read_from_bitstream(uint32_t nitems, uint32_t bits_per_item, uint32_t *unpacked_buffer, const uint32_t *packed_buffer, uint32_t *offset);

//write nitems number of elements, each with bits_per_item bits, from unpacked buffer
// to packed buffer at offset. offset is set to the next aligned uint32_t
void write_to_bitstream(uint32_t nitems, uint32_t bits_per_item, const uint32_t *unpacked_buffer, uint32_t *packed_buffer, uint32_t *offset);

#endif
