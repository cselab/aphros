/**
 * @file checkconsistency.c
 * @authors Simone Bnà, Sandro Manservisi, Ruben Scardovelli,
 *          Philip Yecko and Stephane Zaleski
 * @date  12 November 2015
 * @brief It contains two functions to check the consistency
 *        with a minimum on a cell side and on a cell face;
 *        this greatly reduces the number of calls to the
 *        functions that compute a minimum.
 */

#include "vofi_stddecl.h"

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * check consistency with a minimum in a cell side                            *
 * INPUT: pointer to the implicit function, function value at the two         *
 * endpoints fe, starting point x0, direction sidedir, grid spacing h0        *
 * OUTPUT: f_iat: either the sign to have a positive function at endpoints    *
 * with a minimum inside, or zero if there is no minimum inside (or a zero    *
 * at both endpoints)                                                         *
 * -------------------------------------------------------------------------- */

int vofi_check_side_consistency(
    integrand impl_func, vofi_creal fe[], vofi_creal x0[], vofi_creal sidedir[],
    vofi_creal h0) {
  int i, f_iat, ftmp;
  vofi_real xs[NDIM], dh, f0, f1, fs, ft;

  fs = fe[0] + fe[1];
  if (fs > 0.)
    f_iat = 1;
  else if (fs < 0.)
    f_iat = -1;
  else /* interface through both endpoints */
    f_iat = 0;

  /* for a minimum, the function should decrease from MIN(|fe|)
                                                         towards the interior */
  if (f_iat != 0) {
    dh = MAX(EPS_M * h0, EPS_R);
    f0 = fabs(fe[0]);
    f1 = fabs(fe[1]);
    if (f0 <= f1)
      ft = f0;
    else {
      ft = f1;
      dh = h0 - dh;
    }
    for (i = 0; i < NDIM; i++)
      xs[i] = x0[i] + dh * sidedir[i];
    fs = f_iat * impl_func(xs);
    ftmp = f_iat;
    if (fs >= ft) f_iat = 0;

    /* extra check for very high curvature and multiple intersections */
    if (!f_iat) {
      for (i = 0; i < NDIM; i++)
        xs[i] = x0[i] + 0.5 * h0 * sidedir[i];
      fs = ftmp * impl_func(xs);
      if (fs < ft) f_iat = ftmp;
    }
  }

  return f_iat;
}

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

chk_data vofi_check_face_consistency(
    integrand impl_func, vofi_creal fv[], vofi_creal x0[], vofi_creal sdir[],
    vofi_creal tdir[], vofi_creal h0) {
  int i, iss, ist, iat;
  vofi_real x1[NDIM], xs[NDIM], xt[NDIM], fl[NVER], f0, fs, ft, dh0;
  chk_data ivga;

  f0 = 0.;
  ivga.ivs = ivga.ivt = ivga.igs = ivga.igt = 0;
  for (i = 0; i < NVER; i++)
    f0 += fv[i];

  if (f0 > 0.)
    ivga.iat = 1;
  else if (f0 < 0.)
    ivga.iat = -1;
  else /* interface through four end points */
    ivga.iat = 0;

  /* for a minimum, the function should decrease from MIN(|fv|)
                                                         towards the interior */
  if (ivga.iat != 0) {
    dh0 = MAX(EPS_M * h0, EPS_R);
    f0 = fabs(f0);
    for (i = 0; i < NVER; i++)
      fl[i] = fabs(fv[i]);

    if (fl[0] < f0) {
      f0 = fl[0];
      ist = iss = 1;
    }
    if (fl[1] < f0) {
      f0 = fl[1];
      ivga.ivt = 1;
      iss = 1;
      ist = -1;
    }
    if (fl[2] < f0) {
      f0 = fl[2];
      ivga.ivs = 1;
      ivga.ivt = 0;
      iss = -1;
      ist = 1;
    }
    if (fl[3] < f0) {
      f0 = fl[3];
      ivga.ivt = ivga.ivs = 1;
      ist = iss = -1;
    }

    iat = 0;
    for (i = 0; i < NDIM; i++) {
      x1[i] = x0[i] + h0 * (ivga.ivs * sdir[i] + ivga.ivt * tdir[i]);
      xs[i] = x1[i] + dh0 * iss * sdir[i];
      xt[i] = x1[i] + dh0 * ist * tdir[i];
    }

    fs = ivga.iat * impl_func(xs);
    if (fs < f0) {
      iat = ivga.iat;
      ivga.igs = 1;
    }

    ft = ivga.iat * impl_func(xt);
    if (ft < f0) {
      iat = ivga.iat;
      ivga.igt = 1;
    }

    ivga.iat = iat;
  }

  return ivga;
}
