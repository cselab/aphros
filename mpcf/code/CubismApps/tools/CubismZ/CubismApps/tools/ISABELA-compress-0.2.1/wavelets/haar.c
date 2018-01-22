/*
 * Author: Sriram Lakshminarasimhan <slakshm2 AT ncsu.edu>
 */

#include "haar.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

/*
 * Given an array of num_elements double precision values, perform Haar
 * wavelet transformation on the values. No coefficients are thrown away.
 * Note that the transformation occurs in-place, i.e. buffer is updated.
 */
void haar_wavelet_transform_forward_double (double *buffer, uint32_t num_elements)
{ 
    // Should be able evenly divide the buffer into approximate and detail
    // coefficients in each level. Must be a power of two.
    assert ((num_elements & (num_elements - 1)) == 0);

    double *buffer_tmp = (double *) malloc (num_elements * sizeof (double));
    uint32_t w = num_elements;
    uint32_t i = 0;

    const double sqrt_two = sqrt (2.0);

    // One level at a time
    while (w > 0) {
         w = (w >> 1);

         for(i = 0; i < w; i ++) {
            buffer_tmp [i]     = (buffer [2 * i] + buffer [2 *i + 1]) / sqrt_two;
            buffer_tmp [i + w] = (buffer [2 * i] - buffer [2 *i + 1]) / sqrt_two;
         }

         for(i = 0; i < (w << 1) ; i ++) {
            buffer [i] = buffer_tmp [i];
         }
    }    

    free (buffer_tmp);
}

/*
 * Given an array of num_elements of Haar-transformed double precision data,
 * perform inverse transformation to get back the original array.  Note that the
 * inversion occurs in-place, i.e. buffer is updated.
 */
void haar_wavelet_transform_inverse_double (double *buffer, uint32_t num_elements)
{
    assert ((num_elements & (num_elements - 1)) == 0);

    double *buffer_tmp = (double *) malloc (num_elements * sizeof (double));
    uint32_t w = 1;
    uint32_t i = 0;

    const double sqrt_two = sqrt (2.0);

    while (w < num_elements) {

         for(i = 0; i < w; i ++) {
            buffer_tmp [i * 2]     = (buffer [i] + buffer [i + w]) / sqrt_two;
            buffer_tmp [i * 2 + 1] = (buffer [i] - buffer [i + w]) / sqrt_two;
         }

         for(i = 0; i < (w << 1) ; i ++) {
            buffer [i] = buffer_tmp [i];
         }

         w = (w << 1);
    }

    free (buffer_tmp);
}

/*
 * Single precision version of haar_wavelet_transform_forward_double
 */
void haar_wavelet_transform_forward_float (float *buffer, uint32_t num_elements)
{
    assert ((num_elements & (num_elements - 1)) == 0);

    float *buffer_tmp = (float *) malloc (num_elements * sizeof (float));
    uint32_t w = num_elements;
    uint32_t i = 0;

    const float sqrt_two = sqrt (2.0);

    while (w > 0) {
         w = (w >> 1);

         for(i = 0; i < w; i ++) {
            buffer_tmp [i]     = (buffer [2 * i] + buffer [2 *i + 1]) / sqrt_two;
            buffer_tmp [i + w] = (buffer [2 * i] - buffer [2 *i + 1]) / sqrt_two;
         }

         for(i = 0; i < (w << 1) ; i ++) {
             buffer [i] = buffer_tmp [i];
         }
    }    

    free (buffer_tmp);
}

/*
 * Single precision version of haar_wavelet_transform_inverse_double
 */
void haar_wavelet_transform_inverse_float (float *buffer, uint32_t num_elements)
{
    assert ((num_elements & (num_elements - 1)) == 0);

    float *buffer_tmp = (float *) malloc (num_elements * sizeof (float));
    uint32_t w = 1;
    uint32_t i = 0;

    const float sqrt_two = sqrt (2.0);

    while (w < num_elements) {

         for(i = 0; i < w; i ++) {
            buffer_tmp [(i << 1)]     = (buffer [i] + buffer [i + w]) / sqrt_two;
            buffer_tmp [(i << 1) + 1] = (buffer [i] - buffer [i + w]) / sqrt_two;
         }

         for(i = 0; i < (w << 1) ; i ++) {
             buffer [i] = buffer_tmp [i];
         }

         w = (w << 1);
    }

    free (buffer_tmp);
}

/*
int main (int argc, char *argv [])
{
    // if (argc != 2) {
    //     printf ("Usage: %s <filename>\n", argv [0]);
    //     exit (EXIT_FAILURE);
    // }

    int i = 0;
    double arr [1024];
    double arr_copy [1024];
    double inv [1024];

    int n = 1024;

    for (i = 0; i < n; i ++) {
        arr [i] = (rand () * 15.0 / RAND_MAX);
        arr_copy [i] = arr [i];
    }
   
    char compressed_buffer [65536];
    uint32_t compressed_buffer_sz = 0;

    haar_wavelet_transform_forward_double (arr, 1024);
    memcpy (inv, arr, n * sizeof (double));

    haar_wavelet_transform_inverse_double (arr, 1024);

    for (i = 0; i < n; i ++) {
        printf ("Original = %lf Transformed = %lf Inverted = %lf\n", arr_copy [i], inv [i], arr [i]);
    }
    
    return 0;
}
*/
