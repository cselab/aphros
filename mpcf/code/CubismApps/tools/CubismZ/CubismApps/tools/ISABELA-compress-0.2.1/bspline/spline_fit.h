#ifndef __SPLINE_FIT_H
#define __SPLINE_FIT_H

#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <string>
#include <cstdlib>
#include <stdint.h>

// C linkage
#ifdef __cplusplus
extern "C" {
#endif
void compress_bspline_double (double *buffer, 
                                uint32_t num_elements, 
                                double *coefficients, 
                                uint32_t num_coefficients
                              );

void compress_bspline_float (float *buffer, 
                                uint32_t num_elements, 
                                float *coefficients, 
                                uint32_t num_coefficients
                              );

void decompress_bspline_double (double *buffer, 
                                    uint32_t num_elements, 
                                    double *coefficients, 
                                    uint32_t num_coefficients
                                );

void decompress_bspline_float (float *buffer, 
                                uint32_t num_elements, 
                                float *coefficients, 
                                uint32_t num_coefficients
                              );
#ifdef __cplusplus
}
#endif

#endif
