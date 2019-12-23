/**
 * @file vofi.h
 * @authors Simone Bnà, Sandro Manservisi, Ruben Scardovelli,
 *          Philip Yecko and Stephane Zaleski
 * @date  12 November 2015
 * @brief Prototypes for the vofi library.
 *
 * For any further information on the algorithm implemented in the code
 * @see http://www.sciencedirect.com/science/article/pii/S0010465515004087
 * @see http://www.sciencedirect.com/science/article/pii/S0045793014001480
 */

#ifndef VOFI_H
#define VOFI_H

typedef double (*integrand)(const double[]);

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Starting from point x0 get a zero of the implicit function given by
 * the user, using gradient ascent/descent, then compute its absolute value at a
 * distance hb along the normal direction.
 * @param impl_func pointer to the implicit function
 * @param x0 starting point
 * @param h0 grid spacing
 * @param ndim0 space dimension
 * @param ix0 switch for @p x0 (ix0=1: point x0 is given; ix0=0: use the default
 * value for x0)
 * @param fh "characteristic" function value
 * @note C/C++ API
 */
double vofi_Get_fh(integrand, const double[], double, int, int);

/**
 * @brief Driver to compute the volume fraction value in a given cell in two
 * and three dimensions.
 * @param impl_func pointer to the implicit function
 * @param x0 starting point
 * @param h0 grid spacing
 * @param fh characteristic function value
 * @param ndim0 space dimension
 * @param cc volume fraction value
 * @note C/C++ API
 */
double vofi_Get_cc(integrand, const double[], double, double, int);

#ifdef __cplusplus
}
#endif

#endif
