/*
 * Author: Sriram Lakshminarasimhan <slakshm2 AT ncsu.edu>
 */

/*
 * Given the compressed_buffer, get back the original_buffer (window of
 * floating-point values by inverting the wavelet transformed data.
 * original_buffer is updated and it must be allocated space for atleast
 * num_elements before calling this function.
 */
void decompress_wavelets_float (float *original_buffer, 
                                    uint32_t    num_elements, 
                                    uchar       *compressed_buffer, 
                                    uint32_t    *compressed_buffer_size
                                );

/* 
 * Given an input float precision array, Haar wavelet transform is  
 * applied to convert it to a set of wavelet coefficients. The quantization is
 * done by making [n - ncoefficients] lowest absolute value coefficients 0. 
 * The resulting coefficients and the position of those coefficients are saved
 * within the compressed buffer
 */
void compress_wavelets_float (const float *original_buffer,
                                uint32_t     num_elements,
                                uint32_t     num_coefficients,
                                uchar        *compressed_buffer,
                                uint32_t     *compressed_buffer_size
                             );

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
                                );

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
                             );

