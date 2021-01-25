/**
 * @file getcc.c
 * @authors Simone Bnà, Sandro Manservisi, Ruben Scardovelli,
 *          Philip Yecko and Stephane Zaleski
 * @date  12 November 2015
 * @brief Driver to compute the integration limits and the volume fraction
 *        in two and three dimensions.
 */

#include "vofi_stddecl.h"

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * Driver to compute the volume fraction value in a given cell in two and     *
 * three dimensions                                                           *
 * INPUT:  pointer to the implicit function, starting point x0, grid          *
 * spacing h0, characteristic function value fh, space dimension ndim0        *
 * OUTPUT: cc: volume fraction value                                          *
 * -------------------------------------------------------------------------- */

vofi_real vofi_Get_cc(
    integrand impl_func, vofi_creal x0[], vofi_creal h0, vofi_creal fh,
    vofi_cint ndim0) {
  int nsub;
  vofi_real pdir[NDIM], sdir[NDIM], tdir[NDIM];
  vofi_real side[2 * NSEG]; // TODO
  vofi_real cc;
  dir_data icps;

  icps = vofi_get_dirs(impl_func, x0, pdir, sdir, tdir, h0, fh, ndim0);
  if (icps.icc >= 0)
    cc = (vofi_real)icps.icc;
  else {
    nsub = vofi_get_limits(impl_func, x0, side, pdir, sdir, tdir, h0, ndim0);
    if (ndim0 == 2)
      cc = vofi_get_area(impl_func, x0, side, pdir, sdir, h0, nsub, icps.ipt);
    else
      cc = vofi_get_volume(
          impl_func, x0, side, pdir, sdir, tdir, h0, nsub, icps.ipt);
  }

  return cc;
}
