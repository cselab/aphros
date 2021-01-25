/**
 * @file vofi_stddecl.h
 * @authors Simone Bnà, Sandro Manservisi, Ruben Scardovelli,
 *          Philip Yecko and Stephane Zaleski
 * @date  12 November 2015
 * @brief File containing the prototypes of the vofi functions.
 */

#ifndef VOFI_STDDECL_H
#define VOFI_STDDECL_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
#define Extern extern "C"
#else
#define Extern extern
#endif

#define PREFIX(s) s

#if NOUNDERSCORE
#define SUFFIX(s) s
#else
#define SUFFIX(s) s##_
#endif

#define EXPORT(s) EXPORT_(PREFIX(s))
#define EXPORT_(s) SUFFIX(s)

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define SGN0P(a) ((a < 0) ? -1 : 1)
#define Sq(a) ((a) * (a))
#define Sq3(a) (a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
#define Sqd3(a, b)                                                 \
  ((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + \
   (a[2] - b[2]) * (a[2] - b[2]))
#define SHFT4(a, b, c, d) \
  (a) = (b);              \
  (b) = (c);              \
  (c) = (d)
#define CPSF(s, t, f, g) \
  (s) = (t);             \
  (f) = (g)

#define EPS_M 1.5e-07
#define EPS_LOC 1.5e-07
#define EPS_E 5.0e-07
#define EPS_R 1.0e-14
#define EPS_NOT0 1.0e-50
#define NDIM 3
#define NVER 4
#define NEND 2
#define NLSX 3
#define NLSY 3
#define NLSZ 3
#define NSEG 10

typedef double vofi_real;
typedef const double vofi_creal;
typedef const int vofi_cint;
typedef int* const vofi_int_cpt;
typedef double (*integrand)(vofi_creal[]);

/* xval: coordinates of the minimum or where the sign has changed, fval: local
   function value, sval: distance from the starting point, if applicable,
   iat: sign to have f>0, or if = 0 there is no minimum or no sign change */
typedef struct {
  vofi_real xval[NDIM];
  vofi_real fval;
  vofi_real sval;
  int iat;
} min_data;

/* ivs,ivt: indices to locate the vertex in the face, igs,igt: on/off indices
   for the gradient components, iat: sign to have f>0, or if = 0 there is no
   minimum or no sign change */
typedef struct {
  int ivs;
  int ivt;
  int igs;
  int igt;
  int iat;
} chk_data;

/* icc: full/empty/cut cell (1/0/-1); ipt: tentative number of integration
   points; isb: number of subdivisions, not yet implemented (hence: 0/1) */
typedef struct {
  int icc;
  int ipt;
  int isb;
} dir_data;

/* function prototypes */

/* Fortran APIs */
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
 * @return fh "characteristic" function value
 * @note Fortran API
 */
vofi_real EXPORT(vofi_get_fh)(
    integrand, vofi_creal[], vofi_creal*, vofi_cint*, vofi_cint*);

/**
 * @brief Driver to compute the volume fraction value in a given cell in two
 * and three dimensions.
 * @param impl_func pointer to the implicit function
 * @param x0 starting point
 * @param h0 grid spacing
 * @param fh characteristic function value
 * @param ndim0 space dimension
 * @return cc volume fraction value
 * @note Fortran API
 */
vofi_real EXPORT(vofi_get_cc)(
    integrand, vofi_creal[], vofi_creal*, vofi_creal*, vofi_cint*);

/**
 * @brief compute the zero in a given segment of length s0, the zero is strictly
 * bounded, i.e. f(0)*f(s0) < 0.
 * @param impl_func pointer to the implicit function
 * @param fe function value at the endpoints
 * @param x0 grid spacing
 * @param dir direction
 * @param s0 segment length
 * @param f_sign sign attribute
 * @return sz: length of the segment where f is negative
 */
vofi_real vofi_get_segment_zero(
    integrand, vofi_creal[], vofi_creal[], vofi_creal[], vofi_creal, vofi_cint);

/**
 * @brief check consistency with a minimum in a cell side.
 * @param impl_func pointer to the implicit function
 * @param fe function value at the two endpoints
 * @param x0 grid spacing
 * @param x0 starting point
 * @param sidedir direction
 * @param h0 grid spacing
 * @return f_iat either the sign to have a positive function at endpoints
 * with a minimum inside, or zero if there is no minimum inside (or a zero
 * at both endpoints)
 */
int vofi_check_side_consistency(
    integrand, vofi_creal[], vofi_creal[], vofi_creal[], vofi_creal);

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * check consistency with a minimum in a cell face                            *
 * INPUT: pointer to the implicit function, function value at the four        *
 * vertices fv, starting point x0, secondary and tertiary directions sdir and *
 * tdir, grid spacing h0                                                      *
 * OUTPUT: structure ivga: 2 indices to locate the vertex with the minimum    *
 * function value, 2 on/off indices for the gradient components, iat same as  *
 * f_iat in the previous function                                             *
 * -------------------------------------------------------------------------- */

/**
 * @brief check consistency with a minimum in a cell face.
 * @param impl_func pointer to the implicit function
 * @param fv function value at the four vertices
 * @param x0 starting point
 * @param sdir secondary direction
 * @param tdir tertiary direction
 * @param h0 grid spacing
 * @return structure ivga: 2 indices to locate the vertex with the minimum
 * function value, 2 on/off indices for the gradient components, iat same as
 * @return f_iat either the sign to have a positive function at endpoints
 * with a minimum inside, or zero if there is no minimum inside (or a zero
 * at both endpoints)
 */
chk_data vofi_check_face_consistency(
    integrand, vofi_creal[], vofi_creal[], vofi_creal[], vofi_creal[],
    vofi_creal);

/**
 * @brief compute the value of the implicit function on a 3x3(x3) local grid.
 * @param impl_func pointer to the implicit function
 * @param x0 starting point
 * @param pdir primary direction
 * @param sdir secondary direction
 * @param tdir tertiary direction
 * @param h0 grid spacing
 * @param fh characteristic function value
 * @param ndim0 space dimension
 * @return icps: icc: full/empty/cut cell (1/0/-1); ipt: tentative number
 * of integration points; isb: number of subdivisions, not yet implemented
 * (hence: 0/1)
 */
dir_data vofi_get_dirs(
    integrand, vofi_creal[], vofi_real[], vofi_real[], vofi_real[], vofi_creal,
    vofi_creal, vofi_cint);

/**
 * @brief subdivide the side along the secondary/tertiary (2/3) direction to
 * define rectangles/rectangular hexahedra with or without the interface
 * @param impl_func pointer to the implicit function
 * @param x0 starting point
 * @param pdir primary direction
 * @param lim_intg start/end of each subdivision (lim_intg[0] = 0,
 * lim_intg[nsub] = h0)
 * @param sdir secondary direction
 * @param tdir tertiary direction
 * @param h0 grid spacing
 * @param stdir subdivision direction (2/3)
 * @param ndim0 space dimension
 * @return nsub: total number of subdivisions
 */
int vofi_get_limits(
    integrand, vofi_creal[], vofi_real[], vofi_creal[], vofi_creal[],
    vofi_creal[], vofi_creal, vofi_cint);

/**
 * @brief compute the interface intersections, if any, with a given cell side;
 * these are new internal/external limits of integration.
 * @param impl_func pointer to the implicit function
 * @param fe function value at the endpoints
 * @param x0 starting point
 * @param stdir direction
 * @param h0 grid spacing
 * @param nsub updated number of subdivisions
 * @param lim_intg updated start of new subdivisions
 */
void vofi_get_side_intersections(
    integrand, vofi_real[], vofi_creal[], vofi_real[], vofi_creal[], vofi_creal,
    vofi_int_cpt);

/**
 * @brief get the external limits of integration that are inside the face.
 * these are new internal/external limits of integration.
 * @param impl_func pointer to the implicit function
 * @param xfsa structure with point position with negative f value and function
 * sign attribute
 * @param x0 starting point
 * @param sdir secondary direction
 * @param tdir tertiary direction
 * @param h0 grid spacing
 * @param nsub updated number of subdivisions
 * @param lim_intg updated start of new subdivisions
 */
void vofi_get_face_intersections(
    integrand, min_data, vofi_creal[], vofi_real[], vofi_creal[], vofi_creal[],
    vofi_creal, vofi_int_cpt);

/**
 * @brief compute the function minimum in a given segment of length s0, the
 * search is immediately stopped if a sign change is detected
 * @param impl_func pointer to the implicit function
 * @param fe function value f at the endpoints
 * @param x0 starting point
 * @param dir direction
 * @param s0 segment length
 * @param f_sign sign attribute (to have f > 0 at endpoints)
 * @param max_iter max number of iterations
 * @return xfsa: structure with position, function value, distance from x0
 * and attribute (= 1 if a sign change has been detected) of the minimum or
 * of a point with a different function sign
 */
min_data vofi_get_segment_min(
    integrand, vofi_creal[], vofi_creal[], vofi_creal[], vofi_creal, vofi_cint,
    vofi_cint);

/**
 * @brief compute the function minimum in a cell face, the search is immediately
 * stopped if a sign change is detected.
 * @param impl_func pointer to the implicit function
 * @param x0 starting point
 * @param sdir secondary direction
 * @param tdir tertiary direction
 * @param ivga structure with indices
 * @param h0 grid spacing
 * @return xfsa: structure with position, function value, distance from x0
 * and attribute (= 1 if a sign change has been detected) of the minimum or
 * of a point with a different function sign
 */
min_data vofi_get_face_min(
    integrand, vofi_creal[], vofi_creal[], vofi_creal[], chk_data, vofi_creal);

/**
 * @brief subdivide the side along the secondary/tertiary (2/3) direction to
 * define rectangles/rectangular hexahedra with or without the interface
 * @param impl_func pointer to the implicit function
 * @param x0 starting point
 * @param int_lim_intg internal limits of integration
 * @param pdir primary direction
 * @param sdir secondary direction
 * @param h0 grid spacing
 * @param nintsub number of internal subdivisions
 * @param nintpt tentative number of internal integration points
 * @return area: normalized value of the cut area or 2D volume fraction
 */
vofi_real vofi_get_area(
    integrand, vofi_creal[], vofi_creal[], vofi_creal[], vofi_creal[],
    vofi_creal, vofi_cint, vofi_cint);

/**
 * @brief compute the normalized cut volume with a double Gauss-Legendre
 * quadrature
 * @param impl_func pointer to the implicit function
 * @param x0 starting point
 * @param ext_lim_intg external limits of integration
 * @param pdir primary direction
 * @param sdir secondary direction
 * @param tdir tertiary direction
 * @param h0 grid spacing
 * @param nextsub number of external subdivisions
 * @param nintpt tentative number of internal integration points
 * @return vol: normalized value of the cut volume or 3D volume fraction
 */
vofi_real vofi_get_volume(
    integrand, vofi_creal[], vofi_creal[], vofi_creal[], vofi_creal[],
    vofi_creal[], vofi_creal, vofi_cint, vofi_cint);

#endif
