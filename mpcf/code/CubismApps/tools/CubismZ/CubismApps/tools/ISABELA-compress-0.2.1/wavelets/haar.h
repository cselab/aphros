/*
 * Author: Sriram Lakshminarasimhan <slakshm2 AT ncsu.edu>
 */

#ifndef HAAR_H
#define HAAR_H

#include <stdint.h>

/*
 * Given an array of num_elements double precision values, perform Haar
 * wavelet transformation on the values. No coefficients are thrown away.
 * Note that the transformation occurs in-place, i.e. buffer is updated.
 */
void haar_wavelet_transform_forward_double (double *buffer, uint32_t num_elements);

/*
 * Given an array of num_elements of Haar-transformed double precision data,
 * perform inverse transformation to get back the original array.  Note that the
 * inversion occurs in-place, i.e. buffer is updated.
 */
void haar_wavelet_transform_inverse_double (double *buffer, uint32_t num_elements);

/*
 * Single precision version of haar_wavelet_transform_forward_double
 */
void haar_wavelet_transform_forward_float (float *buffer, uint32_t num_elements);

/*
 * Single precision version of haar_wavelet_transform_inverse_double
 */
void haar_wavelet_transform_inverse_float (float *buffer, uint32_t num_elements);

#endif
