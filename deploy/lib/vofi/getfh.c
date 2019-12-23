/**
 * @file getfh.c
 * @authors Simone Bnà, Sandro Manservisi, Ruben Scardovelli,
 *          Philip Yecko and Stephane Zaleski
 * @date  12 November 2015
 * @brief It computes the characteristic function value fh.
 */

#include "vofi_stddecl.h"

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * Starting from point x0 get a zero of the implicit function given by the    *
 * user, using gradient ascent/descent, then compute its absolute value at a  *
 * distance hb along the normal direction                                     *
 * INPUT:  pointer to the implicit function, starting point x0, grid spacing  *
 * h0, space dimension ndim0, switch ix0 for x0 (ix0=1: point x0 is given;    *
 * ix0=0: use the default value for x0)                                       *
 * OUTPUT: fh ("characteristic" function value)                               *
 * -------------------------------------------------------------------------- */

vofi_real vofi_Get_fh(
    integrand impl_func, vofi_creal x0[], vofi_creal h0, vofi_cint ndim0,
    vofi_cint ix0) {
  int i, k, isw;
  vofi_cint kmax = 100; /* max number of iterations    */
  vofi_creal gamma = 0.01; /* min step along the gradient */
  vofi_creal dh = 1.e-5; /* for 1st deriv. with c.f.d.  */
  vofi_real x1[NDIM], x2[NDIM], xn1[NDIM], xn2[NDIM], der[NDIM], fe[NEND];
  vofi_real f1, f2, fh, delta, dd, hb;

  fh = 4. * h0; /* default value of fh */
  isw = 1;
  if (ndim0 == 3)
    hb = 0.355 * h0;
  else if (ndim0 == 2) {
    hb = 0.255 * h0;
    x1[2] = x2[2] = xn1[2] = xn2[2] = der[2] = 0.;
  } else { /* wrong dimensions! */
    fprintf(stderr, "Wrong dimensions: n =%2d! \n", ndim0);
    fh = -1.;
    isw = 0;
  }

  if (isw) {
    if (ix0 != 0) /* starting point */
      for (i = 0; i < ndim0; i++)
        x2[i] = x0[i];
    else
      for (i = 0; i < ndim0; i++)
        x2[i] = 0.5;

    k = 0;
    f2 = impl_func(x2); /* its f value (should not be zero) */
    while (fabs(f2) < EPS_NOT0 && k < kmax) {
      for (i = 0; i < ndim0; i++)
        x2[i] += dh;
      f2 = impl_func(x2);
      k++;
    }
    f1 = f2;
    if (f1 > 0.) isw = -1;

    /* try to get 2 points with opposite sign of f, by moving along the
                                           gradient direction with step delta */
    k = 0;
    while (f1 * f2 >= 0. && k < kmax) {
      for (i = 0; i < ndim0; i++)
        xn2[i] = xn1[i] = x1[i] = x2[i];
      f1 = f2;
      for (i = 0; i < ndim0; i++) {
        xn2[i] += dh;
        xn1[i] -= dh;
        der[i] = 0.5 * (impl_func(xn2) - impl_func(xn1)) / dh;
        xn2[i] = xn1[i] = x1[i];
      }
      /* DEBUG 1 */

      delta = sqrt(Sq3(der));
      if (delta < EPS_M) {
        for (i = 0; i < ndim0; i++)
          der[i] = 1.;
        delta = sqrt(Sq3(der));
      }
      for (i = 0; i < ndim0; i++)
        der[i] = der[i] / delta;
      delta = fabs(f1 / delta);
      delta = MAX(delta, gamma);
      for (i = 0; i < ndim0; i++)
        x2[i] = x1[i] + isw * delta * der[i];
      f2 = impl_func(x2);
      k++;
    }
    /* DEBUG 2 */

    /* if k<kmax (same as f1*f2 < 0), get the zero on the segment */
    if (k < kmax) {
      delta = sqrt(Sqd3(x1, x2) + EPS_NOT0);
      fe[0] = f1;
      fe[1] = f2;
      for (i = 0; i < ndim0; i++)
        der[i] = (x2[i] - x1[i]) / delta;
      dd = vofi_get_segment_zero(impl_func, fe, x1, der, delta, 1);
      for (i = 0; i < ndim0; i++) {
        if (f1 <= f2)
          x1[i] = x1[i] + dd * der[i];
        else
          x1[i] = x2[i] - dd * der[i];
        xn2[i] = xn1[i] = x1[i];
      }

      /* then get the f value at the distance hb from the zero */
      for (i = 0; i < ndim0; i++) {
        xn2[i] += dh;
        xn1[i] -= dh;
        der[i] = 0.5 * (impl_func(xn2) - impl_func(xn1)) / dh;
        xn2[i] = xn1[i] = x1[i];
      }
      dd = sqrt(Sq3(der) + EPS_NOT0);
      if (dd < EPS_M) {
        fprintf(stderr, "WARNING: the zero is almost a critical point:  \n");
        fprintf(
            stderr, "(x,y,z): (%e, %e, %e) |f|,|grad(f)|: %e, %e \n", x1[0],
            x1[1], x1[2], fabs(impl_func(x1)), dd);
      }
      for (i = 0; i < ndim0; i++) {
        xn1[i] = x1[i] + hb * der[i] / dd;
        xn2[i] = x1[i] - hb * der[i] / dd;
      }
      f1 = fabs(impl_func(xn1));
      f2 = fabs(impl_func(xn2));
      fh = MAX(f1, f2);
    } else { /* did not get f1*f2 < 0! */
      fprintf(stderr, "Did not get f1*f2<0, k,kmax: %3d %3d  \n", k, kmax);
      fh = -1.;
    }
  }
  /* DEBUG 3 */

  return fh;
}
