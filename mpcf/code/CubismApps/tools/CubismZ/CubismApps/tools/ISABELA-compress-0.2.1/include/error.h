/*
 * Author: Sriram Lakshminarasimhan <slakshm2 AT ncsu.edu>
 */

#ifndef __ERROR_HEADER__
#define __ERROR_HEADER__

#include <stdio.h>
#include <stdint.h>

#define MAX_LEVEL 5

uint32_t error_encode_multilevel (double original, double approximated, int32_t *total_errors, float error_rate);
double error_decode_multilevel (double approximated, int32_t *total_errors, float error_rate);

#endif
