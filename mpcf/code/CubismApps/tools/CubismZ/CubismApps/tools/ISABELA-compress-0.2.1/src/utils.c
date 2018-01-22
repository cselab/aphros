/*
 * Author: Sriram Lakshminarasimhan <slakshm2 AT ncsu.edu>
 */

#include "bitstream.h"
#include "utils.h"
#include "string.h"
#include <math.h>
#include "gsl/gsl_sort.h"

struct data_rank {
    double data;
    uint32_t rank;
};

int comparator (const void *a, const void *b)
{
   if ((*(struct data_rank *) a).data > (* (struct data_rank *) b).data) return 1;
   return -1;
}

uint32_t calculate_bits_needed (uint32_t n)
{
    int count = 0;
    while (n > 0) {
        count ++;
        n = n >> 1;
    }

    return count;
}

void sort_and_bit_index (void *input, uchar element_byte_size, uint32_t nelements, uint32_t bits_per_index_element, uint32_t *permutation_index, uint32_t *index_offset)
{
    uint32_t element_idx;
    size_t *size_t_index_tmp_buffer = (size_t *) malloc (nelements * sizeof (size_t));

    gsl_sort_index (size_t_index_tmp_buffer, (double *) input, 1, nelements);
    gsl_sort ((double *) input, 1, nelements);

    if (sizeof (size_t) != sizeof (uint32_t)) {
        uint32_t uint_index_tmp_buffer [nelements];
        uint32_t tmp_index_offset = 0;
        for (element_idx = 0; element_idx < nelements; element_idx ++) {
            uint_index_tmp_buffer [element_idx] = size_t_index_tmp_buffer [element_idx];
        }
        write_to_bitstream (nelements, bits_per_index_element, uint_index_tmp_buffer, permutation_index, index_offset);
    } else {
        write_to_bitstream (nelements, bits_per_index_element, (uint32_t *) size_t_index_tmp_buffer, permutation_index, index_offset);
	}

    free (size_t_index_tmp_buffer);
    return;
}

void sort_and_bit_index_float (void *input, uchar element_byte_size, uint32_t nelements, uint32_t bits_per_index_element, uint32_t *permutation_index, uint32_t *index_offset)
{
    uint32_t element_idx;
    size_t *size_t_index_tmp_buffer = (size_t *) malloc (nelements * sizeof (size_t));

    gsl_sort_float_index (size_t_index_tmp_buffer, (float *) input, 1, nelements);
    gsl_sort_float ((float *) input, 1, nelements);

    if (sizeof (size_t) != sizeof (uint32_t)) {
        uint32_t uint_index_tmp_buffer [nelements];
        uint32_t tmp_index_offset = 0;
        for (element_idx = 0; element_idx < nelements; element_idx ++) {
            uint_index_tmp_buffer [element_idx] = size_t_index_tmp_buffer [element_idx];
        }
        write_to_bitstream (nelements, bits_per_index_element, uint_index_tmp_buffer, permutation_index, index_offset);
    } else {
        write_to_bitstream (nelements, bits_per_index_element, (uint32_t *) size_t_index_tmp_buffer, permutation_index, index_offset);
	}

    free (size_t_index_tmp_buffer);
    return;
}
