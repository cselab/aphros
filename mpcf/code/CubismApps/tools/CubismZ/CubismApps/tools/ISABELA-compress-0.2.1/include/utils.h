/*
 * Author: Sriram Lakshminarasimhan <slakshm2 AT ncsu.edu>
 */

#ifndef __UTILS_HEADER__
#define __UTILS_HEADER__

#include <stddef.h>
#include <stdint.h>

#define MALLOC(buf, type, num) {\
    buf = (type *) malloc (num * sizeof (type)); \
    if (buf == NULL) { \
        printf ("[%s:%d] Unable to allocate memory of %d bytes\n", __FUNCTION__, __LINE__, num * sizeof (type)); \
        exit (-1); \
    } \
}

#define FREE(buf) {\
    free (buf); \
    buf = 0; \
}

#define MIN(A, B) {(A < B? A: B)}
#define MAX(A, B) {(A > B? A: B)}

#define SET(type_name, ptr, value) { \
    (* ((type_name *) ptr)) = value; \
}

#define GET(type_name, ptr, variable) { \
    variable = (* ((type_name *) ptr)); \
}

#define GET_ARRAY(type_name, n, dest_ptr, src_ptr) { \
    memcpy (dest_ptr, src_ptr, n * sizeof(type_name)); \
}

#define SET_ARRAY(type_name, n, dest_ptr, src_ptr) { \
    memcpy (dest_ptr, src_ptr, n * sizeof(type_name)); \
}

typedef unsigned char uchar;

void sort_and_bit_index (void *input, uchar element_byte_size, uint32_t nelements, uint32_t bits_per_index_element, uint32_t *permutation_index, uint32_t *index_offset);

void sort_and_bit_index_float (void *input, uchar element_byte_size, uint32_t nelements, uint32_t bits_per_index_element, uint32_t *permutation_index, uint32_t *index_offset);

uint32_t calculate_bits_needed (uint32_t n);

#endif
