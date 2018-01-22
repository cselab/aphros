#include "spline_fit.h"
#include <stdint.h>
#include "BSpline.h"
#include "ReconstructClass.h"

void compress_bspline_double (double *buffer, 
                                uint32_t num_elements, 
                                double *coefficients, 
                                uint32_t num_coefficients
                              )
{
    uint32_t i = 0;

    double *x = new double [num_elements];
    double *y = buffer;

    for(i = 0; i < num_elements; i ++) {
        x [i] = i;
    }

    BSpline<double> spline (x, num_elements, y, 5, 2, num_coefficients - 1);
     
    if(spline.ok()) {
        for (i = 0; i < num_coefficients - 1; i ++) {         
            coefficients [i] = spline.coefficient (i);
        }
        // mean y value stores as last coefficient
        coefficients [num_coefficients - 1] = spline.getMean ();
    }

    delete [] x;
}

void decompress_bspline_double (double *buffer, 
                                    uint32_t num_elements, 
                                    double *coefficients, 
                                    uint32_t num_coefficients
                                )
{
    Reconstruction R (num_elements, 5, 2, num_coefficients - 1, 0, num_elements - 1); 
    double *x = new double [num_elements];
    double *y = buffer;
    int32_t i = 0;

    for (i = 0; i < num_elements; i ++) {
        x [i] = i;
    }
    double mean = coefficients [num_coefficients - 1];

    BSpline<double> spline (x, coefficients, num_elements, 5, 2, num_coefficients - 1, mean); 

    for (i = 0; i < num_elements; i ++) {
        y [i] = spline.evaluate (x [i]);
    }

    delete [] x;

    return ;
}

void compress_bspline_float (float *buffer, 
                                uint32_t num_elements, 
                                float *coefficients, 
                                uint32_t num_coefficients
                              )
{
    uint32_t i = 0;

    float *x = new float [num_elements];
    float *y = buffer;

    for(i = 0; i < num_elements; i ++) {
        x [i] = i;
    }

    BSpline<float> spline (x, num_elements, y, 5, 2, num_coefficients - 1);
     
    if(spline.ok()) {
        for (i = 0; i < num_coefficients - 1; i ++) {         
            coefficients [i] = spline.coefficient (i);
        }
        // mean y value stores as last coefficient
        coefficients [num_coefficients - 1] = spline.getMean ();
    }

    delete [] x;
}

void decompress_bspline_float (float *buffer, 
                                    uint32_t num_elements, 
                                    float *coefficients, 
                                    uint32_t num_coefficients
                                )
{
    Reconstruction R (num_elements, 5, 2, num_coefficients - 1, 0, num_elements - 1); 
    float *x = new float [num_elements];
    float *y = buffer;

    int32_t i = 0;

    for (i = 0; i < num_elements; i ++) {
        x [i] = i;
    }
    float mean = coefficients [num_coefficients - 1];

    BSpline<float> spline (x, coefficients, num_elements, 5, 2, num_coefficients - 1, mean); 

    for (i = 0; i < num_elements; i ++) {
        y [i] = spline.evaluate (x [i]);
    }

    delete [] x;

    return ;
}
