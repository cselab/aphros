/**
 * @file integrate.c
 * @authors Simone Bnà, Sandro Manservisi, Ruben Scardovelli,
 *          Philip Yecko and Stephane Zaleski
 * @date  12 November 2015
 * @brief It contains two functions to call from Fortran the
 *        corresponding C functions.
 */

#include "vofi_GL.h"
#include "vofi_stddecl.h"

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * compute the normalized cut area with a Gauss-Legendre quadrature           *
 * INPUT: pointer to the implicit function, starting point x0, internal       *
 * limits of integration int_lim_intg, primary and secondary directions pdir  *
 * and  sdir, grid spacing h0, number of internal subdivisions nintsub,       *
 * tentative number of internal integration points nintpt                     *
 * OUTPUT: area: normalized value of the cut area or 2D volume fraction       *
 * -------------------------------------------------------------------------- */

double vofi_get_area(
    integrand impl_func, vofi_creal x0[], vofi_creal int_lim_intg[],
    vofi_creal pdir[], vofi_creal sdir[], vofi_creal h0, vofi_cint nintsub,
    vofi_cint nintpt) {
  int i, ns, k, npt, cut_rect;
  vofi_cint true_sign = 1;
  vofi_real x1[NDIM], x20[NDIM], x21[NDIM], fe[NEND];
  vofi_real area, ds, cs, xis, ht, GL_1D;
  vofi_creal *ptinw, *ptinx;

  /* GRAPHICS I */

  area = 0.;
  for (i = 0; i < NDIM; i++)
    x1[i] = x0[i] + pdir[i] * h0;

  /* DEBUG 1 */

  for (ns = 1; ns <= nintsub; ns++) { /* loop over the rectangles */
    ds = int_lim_intg[ns] - int_lim_intg[ns - 1];
    cs = 0.5 * (int_lim_intg[ns] + int_lim_intg[ns - 1]);
    cut_rect = 0;
    for (i = 0; i < NDIM; i++) {
      x20[i] = x0[i] + sdir[i] * cs;
      x21[i] = x1[i] + sdir[i] * cs;
    }
    fe[0] = impl_func(x20);
    fe[1] = impl_func(x21);
    if (fe[0] * fe[1] <= 0.) cut_rect = 1;

    if (!cut_rect) { /* no interface: full/empty rectangle */
      if (fe[0] < 0.0) area += ds * h0;
      /* DEBUG 2 */

    } else { /* cut rectangle: internal numerical integration */
      if (ds < 0.1 * h0)
        npt = 4;
      else if (ds < 0.2 * h0)
        npt = 8;
      else if (ds < 0.4 * h0)
        npt = MIN(nintpt, 12);
      else if (ds < 0.6 * h0)
        npt = MIN(nintpt, 16);
      else
        npt = MIN(nintpt, 20);

      switch (npt) {
        case 4:
          ptinx = csi04;
          ptinw = wgt04;
          break;
        case 8:
          ptinx = csi08;
          ptinw = wgt08;
          break;
        case 12:
          ptinx = csi12;
          ptinw = wgt12;
          break;
        case 16:
          ptinx = csi16;
          ptinw = wgt16;
          break;
        default:
          ptinx = csi20;
          ptinw = wgt20;
          break;
      }

      GL_1D = 0.;
      /* DEBUG 3 */

      for (k = 0; k < npt; k++) {
        xis = cs + 0.5 * ds * (*ptinx);
        for (i = 0; i < NDIM; i++) {
          x20[i] = x0[i] + sdir[i] * xis;
          x21[i] = x1[i] + sdir[i] * xis;
        }
        fe[0] = impl_func(x20);
        fe[1] = impl_func(x21);
        if (fe[0] * fe[1] < 0.)
          ht = vofi_get_segment_zero(impl_func, fe, x20, pdir, h0, true_sign);
        else { /* weird situation with multiple zeroes */
          if (fe[0] + fe[1] < 0.)
            ht = h0;
          else
            ht = 0.;
        }
        /* DEBUG 4 */

        GL_1D += (*ptinw) * ht;
        ptinx++;
        ptinw++;
        /* GRAPHICS II */
      }
      /* GRAPHICS III */

      area += 0.5 * ds * GL_1D;
    }
  }

  area = area / (h0 * h0); /* normalized area value */

  return area;
}

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * compute the normalized cut volume with a double Gauss-Legendre quadrature  *
 * INPUT: pointer to the implicit function, starting point x0, external       *
 * limits of integration ext_lim_intg, primary, secondary and tertiary        *
 * directions pdir, sdir and tdir, grid spacing h0, number of external        *
 * subdivisions nextsub, tentative number of internal integration points      *
 * nintpt                                                                     *
 * OUTPUT: vol: normalized value of the cut volume or 3D volume fraction      *
 * -------------------------------------------------------------------------- */

double vofi_get_volume(
    integrand impl_func, vofi_creal x0[], vofi_creal ext_lim_intg[],
    vofi_creal pdir[], vofi_creal sdir[], vofi_creal tdir[], vofi_creal h0,
    vofi_cint nextsub, vofi_cint nintpt) {
  int i, ns, k, nexpt, cut_hexa, f_iat, nintsub;
  vofi_cint stdir = 2, max_iter = 50;
  vofi_real x1[NDIM], x2[NDIM], x3[NDIM], fe[NEND], int_lim_intg[NSEG];
  vofi_real vol, ds, cs, xis, f1, f2, area_n, GL_1D;
  vofi_creal *ptexw, *ptexx;
  min_data xfsa;

  vol = 0.;

  /* DEBUG 1 */

  for (ns = 1; ns <= nextsub; ns++) { /* loop over the rectangular hexahedra */
    ds = ext_lim_intg[ns] - ext_lim_intg[ns - 1];
    cs = 0.5 * (ext_lim_intg[ns] + ext_lim_intg[ns - 1]);
    cut_hexa = 0;
    for (i = 0; i < NDIM; i++) {
      x1[i] = x0[i] + tdir[i] * cs;
      x2[i] = x1[i] + pdir[i] * h0;
    }
    f1 = impl_func(x1);
    f2 = impl_func(x2);
    if (f1 * f2 <= 0.) cut_hexa = 1;
    if (!cut_hexa) { /* check lower side along secondary direction */
      fe[0] = f1;
      for (i = 0; i < NDIM; i++)
        x3[i] = x1[i] + sdir[i] * h0;
      fe[1] = impl_func(x3);
      if (fe[0] * fe[1] <= 0.)
        cut_hexa = 1;
      else {
        f_iat = vofi_check_side_consistency(impl_func, fe, x1, sdir, h0);
        if (f_iat != 0) {
          xfsa = vofi_get_segment_min(
              impl_func, fe, x1, sdir, h0, f_iat, max_iter);
          cut_hexa = xfsa.iat;
        }
      }
    }
    if (!cut_hexa) { /* check upper side along secondary direction */
      fe[0] = f2;
      for (i = 0; i < NDIM; i++)
        x3[i] = x2[i] + sdir[i] * h0;
      fe[1] = impl_func(x3);
      if (fe[0] * fe[1] <= 0.)
        cut_hexa = 1;
      else {
        f_iat = vofi_check_side_consistency(impl_func, fe, x2, sdir, h0);
        if (f_iat != 0) {
          xfsa = vofi_get_segment_min(
              impl_func, fe, x2, sdir, h0, f_iat, max_iter);
          cut_hexa = xfsa.iat;
        }
      }
    }

    if (!cut_hexa) { /* no interface: full/empty hexahedron */
      if (f1 < 0.) vol += ds;
      /* DEBUG 2 */

    } else { /* cut hexahedron: external numerical integration */
      if (ds < 0.1 * h0) {
        nexpt = 8;
        ptexx = csi08;
        ptexw = wgt08;
      } else if (ds < 0.3 * h0) {
        nexpt = 12;
        ptexx = csi12;
        ptexw = wgt12;
      } else if (ds < 0.5 * h0) {
        nexpt = 16;
        ptexx = csi16;
        ptexw = wgt16;
      } else {
        nexpt = 20;
        ptexx = csi20;
        ptexw = wgt20;
      }
      GL_1D = 0.;
      /* DEBUG 3 */

      for (k = 0; k < nexpt; k++) {
        xis = cs + 0.5 * ds * (*ptexx);
        for (i = 0; i < NDIM; i++)
          x1[i] = x0[i] + tdir[i] * xis;
        nintsub = vofi_get_limits(
            impl_func, x1, int_lim_intg, pdir, sdir, tdir, h0, stdir);
        area_n = vofi_get_area(
            impl_func, x1, int_lim_intg, pdir, sdir, h0, nintsub, nintpt);
        /* DEBUG 4 */

        GL_1D += (*ptexw) * area_n;
        ptexx++;
        ptexw++;
      }
      vol += 0.5 * ds * GL_1D;
    }
  }

  vol = vol / h0; /* normalized volume value */

  return vol;
}
