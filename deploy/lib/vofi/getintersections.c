/**
 * @file getintersections.c
 * @authors Simone Bnà, Sandro Manservisi, Ruben Scardovelli,
 *          Philip Yecko and Stephane Zaleski
 * @date  12 November 2015
 * @brief It contains two functions to compute the interface
 *        intersection(s) with a cell side and inside a face,
 *        these are internal/external limits of integration.
 */

#include "vofi_stddecl.h"

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * compute the interface intersections, if any, with a given cell side;       *
 * these are new internal/external limits of integration                      *
 * INPUT: pointer to the implicit function, function value at the endpoints   *
 * fe, starting point x0, direction stdir, grid spacing h0                    *
 * OUTPUT: nsub: updated number of subdivisions; array lim_intg: updated      *
 * start of new subdivisions                                                  *
 * -------------------------------------------------------------------------- */

void vofi_get_side_intersections(
    integrand impl_func, vofi_real fe[], vofi_creal x0[], vofi_real lim_intg[],
    vofi_creal stdir[], vofi_creal h0, vofi_int_cpt nsub) {
  int f_iat;
  vofi_cint true_sign = 1, max_iter = 50;
  vofi_real dh0, fh0, ss;
  min_data xfsa;

  if (fe[0] * fe[1] < 0.0) {
    dh0 = vofi_get_segment_zero(impl_func, fe, x0, stdir, h0, true_sign);
    if (fe[0] > 0.0) dh0 = h0 - dh0;
    lim_intg[*nsub] = dh0;
    (*nsub)++;
  } else {
    f_iat = vofi_check_side_consistency(impl_func, fe, x0, stdir, h0);
    if (f_iat != 0) {
      xfsa =
          vofi_get_segment_min(impl_func, fe, x0, stdir, h0, f_iat, max_iter);
      if (xfsa.iat != 0) {
        fh0 = fe[1];
        fe[1] = xfsa.fval;
        dh0 = vofi_get_segment_zero(
            impl_func, fe, x0, stdir, xfsa.sval, true_sign);
        if (fe[0] > 0.0 || fe[1] < 0.0) dh0 = xfsa.sval - dh0;
        lim_intg[*nsub] = dh0;
        (*nsub)++;
        ss = h0 - xfsa.sval;
        fe[0] = fe[1];
        fe[1] = fh0;
        dh0 = vofi_get_segment_zero(
            impl_func, fe, xfsa.xval, stdir, ss, true_sign);
        if (fe[0] > 0.0 || fe[1] < 0.0) dh0 = ss - dh0;
        lim_intg[*nsub] = xfsa.sval + dh0;
        (*nsub)++;
      }
    }
  }

  return;
}

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * get the external limits of integration that are inside the face            *
 * INPUT: pointer to the implicit function, structure with point position     *
 * with negative f value and function sign attribute, starting point x0,      *
 * secondary and tertiary directions sdir and tdir, grid spacing h0           *
 * OUTPUT: nsub: updated number of subdivisions; array lim_intg: updated      *
 * start of new subdivisions                                                  *
 * -------------------------------------------------------------------------- */

void vofi_get_face_intersections(
    integrand impl_func, min_data xfsa, vofi_creal x0[], vofi_real lim_intg[],
    vofi_creal sdir[], vofi_creal tdir[], vofi_creal h0, vofi_int_cpt nsub) {
  int i, k, iter, js, jt, not_conv, ipt, ist, f_iat;
  vofi_cint max_iter = 50;
  vofi_real pt0[NDIM], pt1[NDIM], pt2[NDIM], ptt[NDIM], mp0[NDIM], mp1[NDIM];
  vofi_real ss[NDIM], exdir[NDIM], indir[NDIM], fe[NEND];
  vofi_real ss0, ds0, fpt0, sss, sst, ssx, ssy, tol2, normdir, d1, d2, a1, a2;
  vofi_creal tol = EPS_M;

  /* GRAPHICS I */
  tol2 = 2. * tol;
  for (i = 0; i < NDIM; i++)
    pt0[i] = xfsa.xval[i];
  fe[0] = xfsa.fval;
  f_iat = xfsa.iat;

  /* initialize internal direction and get js and jt indices */
  for (i = 0; i < NDIM; i++) {
    indir[i] = sdir[i];
    if (sdir[i] > 0.5) js = i;
    if (tdir[i] > 0.5) jt = i;
  }

  for (i = 0; i < NDIM; i++)
    pt2[i] = pt1[i] = pt0[i];

  /* get zero or boundary point pt2 along secondary direction with ss -> h0 */
  ss0 = x0[js] + h0 - pt0[js];
  pt2[js] = x0[js] + h0;
  fe[1] = f_iat * impl_func(pt2);
  if (fe[1] > 0.) {
    ds0 = vofi_get_segment_zero(impl_func, fe, pt0, indir, ss0, f_iat);
    pt2[js] = pt0[js] + ds0;
  }
  /* DEBUG 1 */

  /* get zero or boundary point pt1 along secondary direction with ss -> 0 */
  ss0 = pt0[js] - x0[js];
  pt1[js] = x0[js];
  indir[js] = -1.;
  fe[1] = f_iat * impl_func(pt1);
  if (fe[1] > 0.) {
    ds0 = vofi_get_segment_zero(impl_func, fe, pt0, indir, ss0, f_iat);
    pt1[js] = pt0[js] - ds0;
  }
  /* DEBUG 2 */

  for (i = 0; i < NDIM; i++)
    pt0[i] = 0.5 * (pt1[i] + pt2[i]); /* starting midpoint pt0 */
  fpt0 = f_iat * impl_func(pt0);
  ss0 = pt2[js] - pt1[js];
  /* DEBUG 3 */

  /* now get the two external limits along the tertiary direction */
  for (k = -1; k <= 1; k = k + 2) {
    /* DEBUG 4 */

    iter = 0;
    not_conv = 1;
    for (i = 0; i < NDIM; i++) /* initialize external direction */
      exdir[i] = tdir[i];
    exdir[jt] = k;
    if (k < 0)
      sst = pt0[jt] - x0[jt];
    else
      sst = x0[jt] + h0 - pt0[jt];
    for (i = 0; i < NDIM; i++) {
      mp1[i] = pt0[i];
      pt1[i] = mp1[i] + sst * exdir[i];
    }
    fe[0] = fpt0;
    fe[1] = f_iat * impl_func(pt1);
    sss = ss0;
    while (not_conv && iter < max_iter) { /* iterative loop for the limit */
      /* DEBUG 5 */

      if (fe[1] > 0.) {
        ds0 = vofi_get_segment_zero(impl_func, fe, mp1, exdir, sst, f_iat);
        sst = ds0;
      }
      /* DEBUG 6 */

      for (i = 0; i < NDIM; i++) {
        pt1[i] = mp1[i] + sst * exdir[i]; /* zero along the secant line */
        mp0[i] = mp1[i];
        mp1[i] = pt2[i] = ptt[i] = pt1[i];
      }
      /* DEBUG 7 */

      /* try to get other zero along the secondary direction */
      ipt = ist = 0;
      ptt[js] += tol;
      fe[0] = f_iat * impl_func(ptt);
      if (fe[0] < 0.) {
        ipt = 1;
        ssx = x0[js] + h0 - ptt[js];
        indir[js] = 1.;
      } else {
        ptt[js] -= tol2;
        fe[0] = f_iat * impl_func(ptt);
        if (fe[0] < 0.) {
          ipt = 1;
          ssx = ptt[js] - x0[js];
          indir[js] = -1.;
        }
      }
      if (ipt) {
        sss = MIN(
            1.2 * sss, ssx); /* get the segment length along secondary dir */
        for (i = 0; i < NDIM; i++)
          pt2[i] = ptt[i] + sss * indir[i];
        fe[1] = f_iat * impl_func(pt2);
        while (fe[1] < 0. && ist < 3 && sss < ssx) {
          sss = MIN(3. * sss, ssx);
          if (ist == 2) sss = ssx;
          for (i = 0; i < NDIM; i++)
            pt2[i] = ptt[i] + sss * indir[i];
          fe[1] = f_iat * impl_func(pt2);
          ist++;
        }
        if (fe[0] * fe[1] < 0.) { /* get other zero along secondary dir */
          ds0 = vofi_get_segment_zero(impl_func, fe, ptt, indir, sss, f_iat);
          for (i = 0; i < NDIM; i++)
            pt2[i] = ptt[i] + ds0 * indir[i];
        }
        /* DEBUG 8 */

        for (i = 0; i < NDIM; i++) /* get midpoint and width */
          mp1[i] = 0.5 * (pt1[i] + pt2[i]);
        fe[0] = f_iat * impl_func(mp1);
        sss = fabs(pt1[js] - pt2[js]);
      }
      /* DEBUG 9 */

      normdir = 0.;
      for (i = 0; i < NDIM; i++) {
        exdir[i] = mp1[i] - mp0[i]; /* secant direction */
        normdir += exdir[i] * exdir[i];
      }
      normdir = sqrt(normdir) + EPS_NOT0;
      for (i = 0; i < NDIM; i++) {
        exdir[i] = exdir[i] / normdir; /* unit secant direction */
        d1 = SGN0P(exdir[i]);
        d2 = fabs(exdir[i]) + EPS_NOT0;
        a1 = (x0[i] - mp1[i]) / (d1 * d2);
        a2 = (x0[i] + h0 - mp1[i]) / (d1 * d2);
        ss[i] = MAX(a1, a2);
      }
      ssy = MIN(ss[0], ss[1]);
      ssy = MIN(ssy, ss[2]);
      sst = MIN(1.2 * sst, ssy);

      if (!ipt || sss < tol2 || sst < EPS_R) { /* convergence criterion */
        not_conv = 0;
        /* DEBUG 10 */

      } else {
        for (i = 0; i < NDIM; i++)
          pt1[i] = mp1[i] + sst * exdir[i];
        fe[1] = f_iat * impl_func(pt1);
        ist = 0; /* get the segment length along secant dir */
        while (fe[1] < 0. && ist < 3 && sst < ssy) {
          sst = MIN(3. * sst, ssy);
          if (ist == 2) sst = ssy;
          for (i = 0; i < NDIM; i++)
            pt1[i] = mp1[i] + sst * exdir[i];
          fe[1] = f_iat * impl_func(pt1);
          ist++;
        }
      }
      iter++;
    }
    lim_intg[*nsub] = mp1[jt] - x0[jt];
    (*nsub)++;
  }

  return;
}
