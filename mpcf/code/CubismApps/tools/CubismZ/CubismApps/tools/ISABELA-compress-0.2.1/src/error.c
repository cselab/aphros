/*
 * Author: Sriram Lakshminarasimhan <slakshm2 AT ncsu.edu>
 */

#include "error.h"
#include <math.h>
#include <float.h>
#include <assert.h>

uint32_t error_encode_multilevel (double original, double reconstructed, int32_t *total_errors, float error_rate)
{
    int k = -1;
    int relative;

    uint32_t precision_multiplier = 100 / error_rate;

    if (fabs (original) <= DBL_EPSILON) {
        total_errors[0] = k + 1;
        return k + 2;
    }

    if (fabs (reconstructed) <= DBL_EPSILON) {
        reconstructed = FLT_MIN; 
    }

    do {

        k++;
        relative = (int) ((original-reconstructed)/original * precision_multiplier);

        if (relative == precision_multiplier) {
            relative = precision_multiplier - 1;
        }

        total_errors [k + 1] = relative;
        reconstructed = reconstructed / (1 - 1.0 * relative / precision_multiplier);

    } while (relative != 0) ;

    total_errors[0] = k + 1;

    return k + 2;
}

double error_decode_multilevel (double reconstructed, int32_t *total_errors, float error_rate)
{
    uint32_t precision_multiplier = 100 / error_rate;

    if (total_errors [0] == 0) {
        return 0;
    }

    if (fabs (reconstructed) <= DBL_EPSILON) {
        reconstructed = FLT_MIN; 
    }

    for (int k = 0; k < total_errors [0]; k ++) {
        reconstructed = reconstructed / (1 - 1.0 * total_errors [k + 1] / precision_multiplier);
    }

    return reconstructed;
}

int test_encode_decode ()
{
    double original[] = {-5.3, 2.7, 8.9, -11.23, 7.85};
    double approximate[] = {-5.13, 2.76, 18.9, -11.123, 7.99};
    double reconstructed [] = {-5.13, 2.76, 18.9, -11.123, 7.99};

    int n = sizeof (original) / sizeof (double);
    int total_errors [1024];

    int i;
    int j;
    int error_sz = 0;
    int tmp_buffer_offset = 0;

    for (i = 0; i < n; i ++) {
        int level = error_encode_multilevel (original [i],
                                                    approximate [i],
                                                    total_errors + error_sz,
                                                    1
                                                    );
        printf ("Deflate level: %d\n", level);
        error_sz += level;
    }

    for (i = 0; i < error_sz; i ++) {
        printf ("%d ", total_errors [i]);
    }

    printf ("\n");

    error_sz = 0;
    for (i = 0; i < n; i ++) {
        reconstructed [i] = error_decode_multilevel (approximate [i], (int *) total_errors + tmp_buffer_offset, 1);
        printf ("Inflate level: %d\n", ((int *) total_errors) [tmp_buffer_offset]);
        tmp_buffer_offset += ((int *) total_errors) [tmp_buffer_offset] + 1;

        printf ("Original: %lf, Reconstructed: %lf, Error: %lf\n", original [i], reconstructed [i], (original [i] - reconstructed [i]) / original [i] * 100);
    }

    return 1;
}
