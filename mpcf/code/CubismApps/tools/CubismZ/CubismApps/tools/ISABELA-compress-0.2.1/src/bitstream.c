/*
 * Author: John Jenkins <jpjenki2 AT ncsu.edu>
 */
#include "bitstream.h"
#include <stdio.h>

#define BITS_PER_INT (8*sizeof(uint32_t))

//retrieve integer values from bitfield buffer, put into result,
// return the offset of the next int in the buffer stream
void read_from_bitstream(uint32_t nitems, uint32_t min_bits_per_item,
        uint32_t *unpacked_buffer, const uint32_t *packed_buffer, uint32_t *offset) {

    //construct bit_mask
    uint32_t mask = 0;
    uint32_t i;
    for (i = 0; i < min_bits_per_item; i++) {
        mask |= 1 << i;
    }   

    int bit_offset = BITS_PER_INT;
    int diff;
    uint32_t value = packed_buffer[*offset];

    for (i = 0; i < nitems; i++) {
        diff = bit_offset - min_bits_per_item;
        if (diff > 0) {
            unpacked_buffer[i] = (value >> diff) & mask;
            bit_offset -= min_bits_per_item;
        }   
        //perfectly aligned
        else if (diff == 0) {
            unpacked_buffer[i] = value & mask;
            if (i < nitems - 1) {
                value = packed_buffer[++(*offset)];
                bit_offset = BITS_PER_INT;
            }   
        }   
        //partly not in this int
        else if (diff < 0) {
            unpacked_buffer[i] = (value << (-diff)) & mask;
            value = packed_buffer[++(*offset)];
            unpacked_buffer[i] |= (value >> (BITS_PER_INT + diff)) & mask;
            bit_offset = BITS_PER_INT + diff;
        }
    }

    if (bit_offset != BITS_PER_INT) {
        (*offset)++;
    }

}

//WARNING :  All the value should fits certain bits
// Otherwise will cause problem
void write_to_bitstream(uint32_t total_nums, uint32_t min_bits_per_item, const uint32_t *unpacked, uint32_t *packed, uint32_t *offset) 
{
    uint32_t bit_offset = BITS_PER_INT; // use 32,16,or 8 bits as packed bit window
    int diff;
    uint32_t i;
    packed[*offset] = 0;
    for (i = 0; i < total_nums; i++) {
        diff = bit_offset - min_bits_per_item;
        if (diff > 0) {
            packed[*offset] |= (unpacked[i] << diff);
            bit_offset -= min_bits_per_item;
        } else if (diff == 0) {
            packed[(*offset)++] |= unpacked[i];
            packed[(*offset)] = 0;
            bit_offset = BITS_PER_INT;
        } else {
            packed[(*offset)++] |= (unpacked[i] >> (-diff));
            packed[*offset] = (unpacked[i] << (BITS_PER_INT + diff));
            bit_offset = BITS_PER_INT + diff;
        }   
    }   
    if (bit_offset != BITS_PER_INT) {
        (*offset)++;
    }   
}
