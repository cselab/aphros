/**
 * @file getdirs.c
 * @authors Simone Bnà, Sandro Manservisi, Ruben Scardovelli,
 *          Philip Yecko and Stephane Zaleski
 * @date  12 November 2015
 * @brief It checks if the cells is either full or empty, if not
 *        it determines the main, second and third coordinate directions.
 */

#include "vofi_stddecl.h"

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * compute the value of the implicit function on a 3x3(x3) local grid, then   *
 * a) if all f values have the same sign and |f| > fh, then if f > 0 nc = 0,  *
 *    else nc = 1;                                                            *
 * b) else compute gradient components and their average value, order coord.  *
 *    directions, compute tentative number of integration points along the    *
 *    secondary direction                                                     *
 * INPUT: pointer to the implicit function, starting point x0, grid spacing   *
 * h0, characteristic function value fh, space dimension ndim0                *
 * OUTPUT: pdir, sdir, tdir: primary, secondary, tertiary coord. directions;  *
 * structure icps: icc: full/empty/cut cell (1/0/-1); ipt: tentative number   *
 * of integration points; isb: number of subdivisions, not yet implemented    *
 * (hence: 0/1)                                                               *
 * -------------------------------------------------------------------------- */

dir_data vofi_get_dirs(
    integrand impl_func, vofi_creal x0[], vofi_real pdir[], vofi_real sdir[],
    vofi_real tdir[], vofi_creal h0, vofi_creal fh, vofi_cint ndim0) {
  int i, j, k, m, n, np1, np0, nmax, kmax, jt, js, jp, npt_with_grad;
  int cpos[NDIM], cneg[NDIM];
  vofi_creal dh = 1.e-5; /* for 1st deriv. with c.f.d. */
  vofi_creal hh = 0.5 * h0;
  vofi_real df0[NLSZ][NLSX][NLSY][NDIM], f0[NLSZ][NLSX][NLSY];
  vofi_real x1[NDIM], x2[NDIM], xp[NDIM], xm[NDIM], gradf_ave[NDIM];
  vofi_real f1, maxomega, minomega, delomega, tmp, denom;
  dir_data icps;

  /* data initialization */
  minomega = 1.;
  maxomega = -1.;
  icps.icc = -1;
  icps.ipt = 0;
  icps.isb = 1;
  jp = js = jt = npt_with_grad = 0;
  np1 = np0 = 0;
  if (ndim0 == 2) {
    nmax = 9;
    kmax = 0;
    x2[2] = 0.;
  } else {
    nmax = 27;
    kmax = NLSZ - 1;
  }
  for (i = 0; i < ndim0; i++)
    x2[i] = x0[i];

  for (k = 0; k <= kmax; k++) /* get f values on local subgrid */
    for (i = 0; i < NLSX; i++)
      for (j = 0; j < NLSY; j++) {
        x1[0] = x2[0] + i * hh;
        x1[1] = x2[1] + j * hh;
        x1[2] = x2[2] + k * hh;
        f1 = impl_func(x1);
        f0[k][i][j] = f1;
        if (fabs(f1) > fh) {
          if (f1 < 0.)
            np1++;
          else
            np0++;
        }
      }

  /* a): check if the cell is full with np1, or empty with np0 */
  if (np1 == nmax) {
    icps.icc = 1;
    icps.isb = 0;
  } else if (np0 == nmax) {
    icps.icc = 0;
    icps.isb = 0;
  }
  /* b): if not compute the average gradient near the interface (when |f|<fh) */
  if (icps.isb == 1) {
    for (i = 0; i < NDIM; i++) {
      gradf_ave[i] = pdir[i] = sdir[i] = tdir[i] = 0.;
      cpos[i] = cneg[i] = 0;
    }

    for (k = 0; k <= kmax; k++)
      for (i = 0; i < NLSX; i++)
        for (j = 0; j < NLSY; j++) {
          if (fabs(f0[k][i][j]) <= fh) {
            df0[k][i][j][2] = 0.;
            x1[0] = x2[0] + i * hh;
            x1[1] = x2[1] + j * hh;
            x1[2] = x2[2] + k * hh;
            for (m = 0; m < NDIM; m++)
              xp[m] = xm[m] = x1[m];
            for (n = 0; n < ndim0; n++) {
              xp[n] += dh;
              xm[n] -= dh;
              df0[k][i][j][n] = 0.5 * (impl_func(xp) - impl_func(xm)) / dh;
              gradf_ave[n] += df0[k][i][j][n];
              xp[n] = xm[n] = x1[n];
              if (df0[k][i][j][n] > 0.)
                cpos[n] = 1;
              else if (df0[k][i][j][n] < 0.)
                cneg[n] = 1;
            }
            npt_with_grad++;
          }
        }

    for (n = 0; n < NDIM; n++)
      gradf_ave[n] = fabs(gradf_ave[n]);

    /* get main, second, and third directions */
    jp = 0;
    js = 1;
    jt = 2;
    if (gradf_ave[1] > gradf_ave[0]) {
      jp = 1;
      js = 0;
    }
    if (gradf_ave[2] > gradf_ave[jp]) {
      jt = js;
      js = jp;
      jp = 2;
    } else if (gradf_ave[2] > gradf_ave[js]) {
      jt = js;
      js = 2;
    }
    pdir[jp] = sdir[js] = tdir[jt] = 1.;

    /* check sign consistency in cut cells with high curvature */
    if (cpos[jp] * cneg[jp] == 1 && cpos[js] * cneg[js] == 0) {
      pdir[jp] = sdir[js] = 0.;
      k = jp;
      jp = js;
      js = k;
      pdir[jp] = sdir[js] = 1.;
    }

    /* get tentative number of integration points along the secondary
                                direction (ns=1: no level of subdivision yet) */
    for (k = 0; k <= kmax; k++)
      for (i = 0; i < NLSX; i++)
        for (j = 0; j < NLSY; j++) {
          if (fabs(f0[k][i][j]) <= fh) {
            denom = Sq(df0[k][i][j][jp]) + Sq(df0[k][i][j][js]) + EPS_NOT0;
            tmp = df0[k][i][j][js] * fabs(df0[k][i][j][js]);
            if (tmp > maxomega * denom) maxomega = tmp / denom;
            if (tmp < minomega * denom) minomega = tmp / denom;
          }
        }
    delomega = maxomega - minomega;
    if (delomega <= 0.35)
      icps.ipt = 8;
    else if (delomega <= 0.65)
      icps.ipt = 12;
    else if (delomega <= 0.9)
      icps.ipt = 16;
    else
      icps.ipt = 20;
  }

  /* DEBUG 1 */

  return icps;
}
