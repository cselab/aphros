/*
 * Author: Sriram Lakshminarasimhan <slakshm2 AT ncsu.edu>
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_statistics.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include "bitstream.h"
#include "haar.h"
#include "utils.h"

#include "openssl/md5.h"

// Serialized view of the the compressed buffer after the compress_wavelet_* ().
// <uint32_t num_total_elements>                 - Total number of elements (with padding) used in transform 
// <uint32_t num_coefficients>                   - Number of coefficients saved
// <void *coefficients>                          - Value of the saved num_coefficients coefficients
// <uint32_t packed_coefficients_position_size>  - Size of the position array of bitpacked num_coefficients coefficients
// <uint32_t *packed_coefficients_position_size> - Bitpacked array of coefficient positions


/*
 * Given the compressed_buffer, get back the original_buffer (window of
 * single-precision values) by inverting the wavelet transformed data.
 * original_buffer is updated and it must be allocated space for atleast
 * num_elements before calling this function.
 * [single-precision version of decompress_wavelets_double]
 */
void decompress_wavelets_float (float *original_buffer, 
                                    uint32_t num_elements, 
                                    uchar *compressed_buffer, 
                                    uint32_t *compressed_buffer_size
                                )
{

    float *wavelet_transform_buffer = 0;
    float *coefficients = 0;

    uint32_t num_total_elements = 0;
    uint32_t num_coefficients = 0;

    uint32_t *coefficients_position = 0;
    uint32_t index_size = 0;
    uint32_t packed_coefficients_position_size = 0;

    uint32_t bits_per_index_element = 0;

    uchar *compressed_buffer_start = compressed_buffer;
    uint32_t i = 0;

    // Get num_pad_elements
    GET (uint32_t, compressed_buffer, num_total_elements)
    compressed_buffer += sizeof (uint32_t);

    // Get num_coefficients
    GET (uint32_t, compressed_buffer, num_coefficients)
    compressed_buffer += sizeof (uint32_t);

    // Get coefficients
    coefficients = (float *) malloc (sizeof (float) * num_coefficients);
    GET_ARRAY (float, num_coefficients, coefficients, compressed_buffer)
    compressed_buffer += num_coefficients * sizeof (float);

    // Get the size of packed coefficients
    GET (uint32_t, compressed_buffer, packed_coefficients_position_size)
    compressed_buffer += sizeof (uint32_t);

    // Unpack coefficients_position_index
    bits_per_index_element = calculate_bits_needed (num_total_elements) - 1;
    coefficients_position = (uint32_t *) malloc (sizeof (uint32_t) * (num_coefficients + 1));
    read_from_bitstream (num_coefficients, bits_per_index_element, coefficients_position, (uint32_t *) (compressed_buffer), &index_size);

    assert (packed_coefficients_position_size == index_size);
    compressed_buffer += packed_coefficients_position_size * sizeof (uint32_t);

    // Place the coefficients in the correct location as indicated by the
    // coefficients_position array
    wavelet_transform_buffer = (float *) malloc (num_total_elements * sizeof (float));
    assert (wavelet_transform_buffer != 0);

    for (i = 0; i < num_total_elements; i ++) {
        wavelet_transform_buffer [i] = 0;
    }

    for (i = 0; i < num_coefficients; i ++) {
        wavelet_transform_buffer [coefficients_position [i]] = coefficients [i];
    }

    haar_wavelet_transform_inverse_float (wavelet_transform_buffer, num_total_elements);

    memcpy (original_buffer, wavelet_transform_buffer, sizeof (float) * num_elements);
        
    *compressed_buffer_size = (compressed_buffer - compressed_buffer_start);

    free (coefficients_position);
    free (wavelet_transform_buffer);
    free (coefficients);

    return;
}

/* 
 * Given an input single-precision array, Haar wavelet transform is  
 * applied to convert it to a set of wavelet coefficients. The quantization is
 * done by making [n - ncoefficients] lowest absolute value coefficients 0. 
 * The resulting coefficients and the position of those coefficients are saved
 * within the compressed buffer.
 * [single-precision floating-point version of compress_wavelets_double]
 */
void compress_wavelets_float (const float *original_buffer,
                                uint32_t     num_elements,
                                uint32_t     num_coefficients,
                                uchar         *compressed_buffer,
                                uint32_t     *compressed_buffer_size
                             )
{
    // If the total number of elements is less than the number of coefficients 
    // try to reduce the number of coefficients used. This should happen only
    // for the last linearized window
    if (num_elements == 1) {
        num_coefficients = 1;
    } else if (num_elements <= num_coefficients) {
        num_coefficients = num_elements / 2;
    }

    // Number of bits to use per index element for marking the position of the
    // saved coefficients.
    uint32_t bits_per_index_element = calculate_bits_needed (num_elements - 1);
    uint32_t num_total_elements = 0;

    if ((num_elements & (num_elements - 1)) == 0) {
        num_total_elements = num_elements;
    } else {
        bits_per_index_element ++; 
        num_total_elements = (1 << bits_per_index_element);
    }

    // If the total number of elements is not a multiple of 2^(WAVELET_NLEVEL),
    // then pad values to the array. The padded value is the maximum value in
    // the array.
    uint32_t num_pad_elements = num_total_elements - num_elements;

    float *wavelet_transform_buffer = (float *) malloc (num_total_elements * sizeof (float));
    float *absolute_coefficients    = (float *) malloc (num_total_elements * sizeof (float));

    uint32_t i = 0;
    uchar *compressed_buffer_start = compressed_buffer;

    size_t *index = (size_t *) malloc (num_total_elements * sizeof (size_t));

    assert (wavelet_transform_buffer != 0);
    assert (absolute_coefficients  != 0);
    assert (index != 0);

    memcpy (wavelet_transform_buffer, original_buffer, num_total_elements * sizeof (float));

    // Pad with the maximum value
    for (i = 0; i < num_pad_elements; i ++) {
        wavelet_transform_buffer [i + num_elements] = wavelet_transform_buffer [num_elements - 1];
    }

    // Perform transformation using HAAR wavelets  
    haar_wavelet_transform_forward_float (wavelet_transform_buffer, num_total_elements);

    // Original buffer has been replaced with wavelet coefficients.
    // Get the coefficients whose absolute value is minimal.
    for (i = 0; i < num_total_elements; i++) {
        absolute_coefficients [i] = fabs (wavelet_transform_buffer [i]);
    }
    gsl_sort_float_index (index, absolute_coefficients, 1, num_total_elements);

    SET (uint32_t, compressed_buffer, num_total_elements)
    compressed_buffer += sizeof (uint32_t);

    SET (uint32_t, compressed_buffer, num_coefficients)
    compressed_buffer += sizeof (uint32_t);

    // Save the top n coefficients
    uint32_t top = 0;
    uint32_t *coefficients_position = (uint32_t *) malloc (sizeof (uint32_t) * num_coefficients);
    uint32_t packed_coefficients_position_size = 0;

    for (i = num_total_elements - num_coefficients; i < num_total_elements; i ++) {
        ((float *) compressed_buffer) [top] = wavelet_transform_buffer [index [i]];
        coefficients_position [top] = index [i];
        top ++;
    }

    // Store the locations of top n coefficients in bit-packed form
    compressed_buffer += num_coefficients * sizeof (float);
    write_to_bitstream (num_coefficients, bits_per_index_element, coefficients_position, (uint32_t *) (compressed_buffer + sizeof (uint32_t)), &packed_coefficients_position_size);

    SET (uint32_t, compressed_buffer, packed_coefficients_position_size);

    compressed_buffer += sizeof (uint32_t);
    compressed_buffer += sizeof (uint32_t) * packed_coefficients_position_size;
    *compressed_buffer_size = (compressed_buffer - compressed_buffer_start);

    // Clear buffers
    free (coefficients_position);
    free (absolute_coefficients);
    free (index);
    free (wavelet_transform_buffer);

    return ;
}

/*
 * Given the compressed_buffer, get back the original_buffer (window of
 * double precision values) by inverting the wavelet transformed data.
 * original_buffer is updated and it must be allocated space for atleast
 * num_elements before calling this function.
 * [double precision version of decompress_wavelets_float]
 */
void decompress_wavelets_double (double *original_buffer, 
                                    uint32_t num_elements, 
                                    uchar *compressed_buffer, 
                                    uint32_t *compressed_buffer_size
                                )
{

    double *wavelet_transform_buffer = 0;
    double *coefficients = 0;

    uint32_t num_total_elements = 0;
    uint32_t num_coefficients = 0;

    uint32_t *coefficients_position = 0;
    uint32_t index_size = 0;
    uint32_t packed_coefficients_position_size = 0;

    uint32_t bits_per_index_element = 0;

    uchar *compressed_buffer_start = compressed_buffer;
    uint32_t i = 0;

    // Get num_pad_elements
    GET (uint32_t, compressed_buffer, num_total_elements)
    compressed_buffer += sizeof (uint32_t);

    // Get num_coefficients
    GET (uint32_t, compressed_buffer, num_coefficients)
    compressed_buffer += sizeof (uint32_t);

    // Get coefficients
    coefficients = (double *) malloc (sizeof (double) * num_coefficients);
    GET_ARRAY (double, num_coefficients, coefficients, compressed_buffer)
    compressed_buffer += num_coefficients * sizeof (double);

    // Get the size of packed coefficients
    GET (uint32_t, compressed_buffer, packed_coefficients_position_size)
    compressed_buffer += sizeof (uint32_t);

    // Unpack coefficients_position_index
    bits_per_index_element = calculate_bits_needed (num_total_elements) - 1;
    coefficients_position = (uint32_t *) malloc (sizeof (uint32_t) * (num_coefficients + 1));
    read_from_bitstream (num_coefficients, bits_per_index_element, coefficients_position, (uint32_t *) (compressed_buffer), &index_size);

    assert (packed_coefficients_position_size == index_size);
    compressed_buffer += packed_coefficients_position_size * sizeof (uint32_t);

    // Place the coefficients in the correct location as indicated by the
    // coefficients_position array
    wavelet_transform_buffer = (double *) malloc (num_total_elements * sizeof (double));
    assert (wavelet_transform_buffer != 0);

    for (i = 0; i < num_total_elements; i ++) {
        wavelet_transform_buffer [i] = 0;
    }

    for (i = 0; i < num_coefficients; i ++) {
        wavelet_transform_buffer [coefficients_position [i]] = coefficients [i];
    }

    haar_wavelet_transform_inverse_double (wavelet_transform_buffer, num_total_elements);

    memcpy (original_buffer, wavelet_transform_buffer, sizeof (double) * num_elements);
        
    *compressed_buffer_size = (compressed_buffer - compressed_buffer_start);

    free (coefficients_position);
    free (wavelet_transform_buffer);
    free (coefficients);

    return;
}

/* 
 * Given an input double precision array, Haar wavelet transform is  
 * applied to convert it to a set of wavelet coefficients. The quantization is
 * done by making [n - ncoefficients] lowest absolute value coefficients 0. 
 * The resulting coefficients and the position of those coefficients are saved
 * within the compressed buffer.
 * [double precision version of compress_wavelets_float]
 */
void compress_wavelets_double (const double *original_buffer,
                                uint32_t     num_elements,
                                uint32_t     num_coefficients,
                                uchar         *compressed_buffer,
                                uint32_t     *compressed_buffer_size
                             )
{
    // If the total number of elements is less than the number of coefficients 
    // try to reduce the number of coefficients used. This should happen only
    // for the last linearized window
    if (num_elements == 1) {
        num_coefficients = 1;
    } else if (num_elements <= num_coefficients) {
        num_coefficients = num_elements / 2;
    }

    // Number of bits to use per index element for marking the position of the
    // saved coefficients.
    uint32_t bits_per_index_element = calculate_bits_needed (num_elements - 1);
    uint32_t num_total_elements = 0;

    if ((num_elements & (num_elements - 1)) == 0) {
        num_total_elements = num_elements;
    } else {
        bits_per_index_element ++; 
        num_total_elements = (1 << bits_per_index_element);
    }

    // If the total number of elements is not a multiple of 2^(WAVELET_NLEVEL),
    // then pad values to the array. The padded value is the maximum value in
    // the array.
    uint32_t num_pad_elements = num_total_elements - num_elements;

    double *wavelet_transform_buffer = (double *) malloc (num_total_elements * sizeof (double));
    double *absolute_coefficients    = (double *) malloc (num_total_elements * sizeof (double));

    uint32_t i = 0;
    uchar *compressed_buffer_start = compressed_buffer;

    size_t *index = (size_t *) malloc (num_total_elements * sizeof (size_t));

    assert (wavelet_transform_buffer != 0);
    assert (absolute_coefficients  != 0);
    assert (index != 0);

    memcpy (wavelet_transform_buffer, original_buffer, num_total_elements * sizeof (double));

    // Pad with the maximum value
    for (i = 0; i < num_pad_elements; i ++) {
        wavelet_transform_buffer [i + num_elements] = wavelet_transform_buffer [num_elements - 1];
    }

    // Perform transformation using HAAR wavelets  
    haar_wavelet_transform_forward_double (wavelet_transform_buffer, num_total_elements);

    // Original buffer has been replaced with wavelet coefficients.
    // Get the coefficients whose absolute value is minimal.
    for (i = 0; i < num_total_elements; i++) {
        absolute_coefficients [i] = fabs (wavelet_transform_buffer [i]);
    }
    gsl_sort_index (index, absolute_coefficients, 1, num_total_elements);

    SET (uint32_t, compressed_buffer, num_total_elements)
    compressed_buffer += sizeof (uint32_t);

    SET (uint32_t, compressed_buffer, num_coefficients)
    compressed_buffer += sizeof (uint32_t);

    // Save the top n coefficients
    uint32_t top = 0;
    uint32_t *coefficients_position = (uint32_t *) malloc (sizeof (uint32_t) * num_coefficients);
    uint32_t packed_coefficients_position_size = 0;

    for (i = num_total_elements - num_coefficients; i < num_total_elements; i ++) {
        ((double *) compressed_buffer) [top] = wavelet_transform_buffer [index [i]];
        coefficients_position [top] = index [i];
        top ++;
    }

    // Store the locations of top n coefficients in bit-packed form
    compressed_buffer += num_coefficients * sizeof (double);
    write_to_bitstream (num_coefficients, bits_per_index_element, coefficients_position, (uint32_t *) (compressed_buffer + sizeof (uint32_t)), &packed_coefficients_position_size);

    SET (uint32_t, compressed_buffer, packed_coefficients_position_size);

    compressed_buffer += sizeof (uint32_t);
    compressed_buffer += sizeof (uint32_t) * packed_coefficients_position_size;
    *compressed_buffer_size = (compressed_buffer - compressed_buffer_start);

    // Clear buffers
    free (coefficients_position);
    free (absolute_coefficients);
    free (index);
    free (wavelet_transform_buffer);

    return ;
}

#ifdef TEST
int test_example (int option)
{
    uint32_t i = 0;
    uint32_t n = 101;

    if (option == 1) {
        double arr [n];
        double arr_copy [n];

        arr [i] = 0;
        for (i = 1; i < n; i ++) {
            arr [i] = arr [i - 1] + (rand () * 15.0 / RAND_MAX);
        }
   
        uchar compressed_buffer [65536];
        uint32_t compressed_buffer_size = 0;

        compress_wavelets_double (arr, n, 30, compressed_buffer, &compressed_buffer_size);
        decompress_wavelets_double (arr_copy, n, compressed_buffer, compressed_buffer_size);

        for (i = 0; i < n; i ++) {
            printf ("%lf %lf\n", arr [i], arr_copy [i]); 
        }
    } else {
        float arr [n];
        float arr_copy [n];

        arr [i] = 0;
        for (i = 1; i < n; i ++) {
            arr [i] = arr [i - 1] + (rand () * 15.0 / RAND_MAX);
        }
   
        uchar compressed_buffer [65536];
        uint32_t compressed_buffer_size = 0;

        compress_wavelets_float (arr, n, 30, compressed_buffer, &compressed_buffer_size);
        decompress_wavelets_float (arr_copy, n, compressed_buffer, compressed_buffer_size);

        for (i = 0; i < n; i ++) {
            printf ("%f %f\n", arr [i], arr_copy [i]); 
        }
    }

    return 0;
}

int main (int argc, char *argv [])
{
    test_example (argc);
}
#endif
