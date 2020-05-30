/**
 * @file getmin.c
 * @authors Simone Bnà, Sandro Manservisi, Ruben Scardovelli,
 *          Philip Yecko and Stephane Zaleski
 * @date  12 November 2015
 * @brief It contains two functions to compute the function minimum
 *        either in a given segment or in a cell face, the search is
 *        stopped if a sign change is detected.
 */

#include "vofi_stddecl.h"

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * compute the function minimum in a given segment of length s0, the search   *
 * is immediately stopped if a sign change is detected                        *
 * METHOD: standard Brent's method plus linear extrapolation for vertices     *
 * INPUT: pointer to the implicit function, function value f at the endpoints *
 * fe, starting point x0, direction dir, segment length s0, sign attribute    *
 * f_sign (to have f > 0 at endpoints), max number of iterations max_iter     *
 * OUTPUT: xfsa: structure with position, function value, distance from x0    *
 * and attribute (= 1 if a sign change has been detected) of the minimum or   *
 * of a point with a different function sign                                  *
 * -------------------------------------------------------------------------- */

min_data vofi_get_segment_min(
    integrand impl_func, vofi_creal fe[], vofi_creal x0[], vofi_creal dir[],
    vofi_creal s0, vofi_cint f_sign, vofi_cint max_iter) {
  int i, j, iter, not_conv, igold, iseca;
  vofi_real xs[NDIM], fs, fu, ft, fv, fa, fb, p, q, r;
  vofi_real sa, sb, ss, su, st, sv, se, sd, sc, sm, fm, sp, fp, sz, fz;
  vofi_real GRIS, tol, t2;
  min_data xfsa;

  /* square of the inverse of the golden ratio */
  GRIS = 0.5 * (3.0 - sqrt(5.0));

  xfsa.iat = 0;
  igold = 1;
  not_conv = 1;

  /* initialize variables to adapt to Brent's routine (f(ss) < f(st) < f(sv)) */
  fa = f_sign * fe[0];
  sa = 0.;
  fb = f_sign * fe[1];
  sb = s0;
  if (fa <= fb) {
    st = 0.;
    sv = s0;
    ft = fa;
    fv = fb;
  } else {
    st = s0;
    sv = 0.;
    ft = fb;
    fv = fa;
  }

  ss = sa + GRIS * (sb - sa);
  for (i = 0; i < NDIM; i++)
    xs[i] = x0[i] + ss * dir[i];
  fs = f_sign * impl_func(xs);
  if (fs > ft) {
    SHFT4(fu, ft, fs, fu);
    SHFT4(su, st, ss, su);
  }
  se = st - sv;
  sd = ss - st;
  if (fs < 0.) /* got the sign change right away! */
    not_conv = 0;

  iter = 0;
  /* DEBUG 1 */

  while (not_conv && iter < max_iter) {
    sc = 0.5 * (sa + sb);
    tol = EPS_M * fabs(ss) + EPS_LOC;
    t2 = 2.0 * tol;
    /* DEBUG 2 */

    /* check if the stopping criterion is satisfied */
    if (fabs(ss - sc) <= t2 - 0.5 * (sb - sa)) not_conv = 0;

    /* else go on and compute the parabola coefficients */
    else {
      iter++;
      r = 0.0;
      q = r;
      p = q;
      if (fabs(se) > tol) {
        r = (ss - st) * (fs - fv);
        q = (ss - sv) * (fs - ft);
        p = (ss - sv) * q - (ss - st) * r;
        q = 2.0 * (q - r);
        if (q > 0.0)
          p = -p;
        else
          q = fabs(q);
        r = se;
        se = sd;
      }
      /* DEBUG 3 */

      /* check if the parabolic interpolation step is OK, but the function
                                  must not be evaluated too close to sa or sb */
      if (fabs(p) < fabs(0.5 * q * r) && p > q * (sa - ss) &&
          p < q * (sb - ss)) {
        sd = p / q;
        su = ss + sd;
        if ((su - sa) < t2 || (sb - su) < t2) {
          if (ss < sc)
            sd = tol;
          else
            sd = -tol;
        }
        /* DEBUG 4 */

      }
      /* otherwise take a golden-section step */
      else {
        if (ss < sc)
          se = sb - ss;
        else
          se = sa - ss;
        sd = GRIS * se;
        igold++;
        /* DEBUG 5 */
      }

      /* the function must not be evaluated too close to ss as well */
      if (fabs(sd) >= tol)
        su = ss + sd;
      else if (sd > 0.0)
        su = ss + tol;
      else
        su = ss - tol;

      /* get new point and f value */
      for (i = 0; i < NDIM; i++)
        xs[i] = x0[i] + su * dir[i];
      fu = f_sign * impl_func(xs);
      if (fu < 0.) /* got the sign change! */
        not_conv = 0;

      /* DEBUG 6 */

      /* update when fu <= fs */
      if (fu <= fs) {
        if (su < ss) {
          CPSF(sb, ss, fb, fs);
        } else {
          CPSF(sa, ss, fa, fs);
        }

        SHFT4(fv, ft, fs, fu);
        SHFT4(sv, st, ss, su);
      }
      /* update when fu > fs */
      else {
        if (su < ss) {
          CPSF(sa, su, fa, fu);
        } else {
          CPSF(sb, su, fb, fu);
        }

        if (fu <= ft || st == ss) {
          CPSF(sv, st, fv, ft);
          CPSF(st, su, ft, fu);
        } else if (fu <= fv || sv == ss || sv == st) {
          CPSF(sv, su, fv, fu);
        }
      }
      /* DEBUG 7 */

      /* every other golden-section step, check a linear extrapolation step */
      if (igold == 2 && not_conv == 1) {
        igold = 0;
        CPSF(sm, ss, fm, fs);
        if (st < sm && fabs(st - sa) > t2) {
          CPSF(sm, st, fm, ft);
        }
        if (sv < sm && fabs(sv - sa) > t2) {
          CPSF(sm, sv, fm, fv);
        }

        CPSF(sp, ss, fp, fs);
        if (st > sp && fabs(st - sb) > t2) {
          CPSF(sp, st, fp, ft);
        }
        if (sv > sp && fabs(sv - sb) > t2) {
          CPSF(sp, sv, fp, fv);
        }

        p = (sa - sm) * (fp * sb - fb * sp) + (sp - sb) * (fm * sa - fa * sm);
        q = (sa - sm) * (fp - fb) + (sp - sb) * (fm - fa);
        if (q < 0.0) {
          p = -p;
          q = -q;
        }

        /* check if the intersection is OK and is a minimum */
        /* DEBUG 8 */

        sm = MIN(sa, sm);
        sp = MAX(sb, sp);
        /* DEBUG 9 */

        if (p > q * sm && p < q * sp) {
          su = p / q;
          for (i = 0; i < NDIM; i++)
            xs[i] = x0[i] + su * dir[i];
          fu = f_sign * impl_func(xs);
          /* DEBUG 10 */

          iseca = 0;
          if (fu < fs) {
            if (fu < 0.) /* got the sign change! */
              not_conv = 0;
            tol = EPS_M * fabs(su) + EPS_LOC;
            for (j = -1; j <= 1; j = j + 2) {
              sz = su + j * tol;
              for (i = 0; i < 3; i++)
                xs[i] = x0[i] + sz * dir[i];
              fz = f_sign * impl_func(xs);
              if (fz > fu) iseca++;
            }
            if (iseca == 2) {
              CPSF(ss, su, fs, fu);
              sb = sz;
              sa = sz - 2. * tol;
              /* DEBUG 11 */
            }
          }
        }
      } /* end extrapolation */
    } /* end iteration */
  } /* end external while */

  /* DEBUG 12 */

  for (i = 0; i < NDIM; i++)
    xfsa.xval[i] = xs[i];
  xfsa.fval = f_sign * fs;
  xfsa.sval = ss;
  if (fs < 0.) xfsa.iat = 1;

  return xfsa;
}

/* -------------------------------------------------------------------------- *
 * DESCRIPTION:                                                               *
 * compute the function minimum in a cell face, the search is immediately     *
 * stopped if a sign change is detected                                       *
 * METHOD: standard preconditioned conjugate gradient method, derivatives     *
 * with centered finite differences, preconditioner is the Hessian main       *
 * diagonal, Polak-Ribiere expression for the beta parameter                  *
 * INPUT: pointer to the implicit function, starting point x0, secondary and  *
 * tertiary directions sdir and tdir, structure with indices ivga, grid       *
 * OUTPUT: xfsa: structure with position, function value, distance from x0    *
 * and attribute (= 1 if a sign change has been detected) of the minimum or   *
 * of a point with a different function sign                                  *
 * -------------------------------------------------------------------------- */

min_data vofi_get_face_min(
    integrand impl_func, vofi_creal x0[], vofi_creal sdir[], vofi_creal tdir[],
    chk_data ivga, vofi_creal h0) {
  int i, not_conv, iter, k, ipt, iss;
  vofi_cint max_iter = 50, max_iter_line = 50;
  vofi_real xs0[NDIM], xs1[NDIM], xs2[NDIM], xt1[NDIM], xt2[NDIM], res[NDIM],
      hes[NDIM];
  vofi_real rs0[NDIM], hs0[NDIM], pcrs[NDIM], nmdr[NDIM], cndr[NDIM], ss[NDIM],
      fe[NEND];
  vofi_real eps2, fp0, fs1, fs2, ft1, ft2, dfs, dft, d2fs, d2ft, mcd, ss0, ss1,
      beta;
  vofi_real del0, delnew, delold, delmid, d1, d2, a1, a2;
  vofi_creal dh = 1.e-04; /* for 1st and 2nd derivatives with c.f.d. */
  min_data xfsa;

  eps2 = EPS_E * EPS_E;

  for (i = 0; i < NDIM; i++) { /* data at starting point */
    xs0[i] = x0[i] + h0 * (ivga.ivs * sdir[i] + ivga.ivt * tdir[i]);
    xs2[i] = xs0[i] + dh * sdir[i];
    xs1[i] = xs0[i] - dh * sdir[i];
    xt2[i] = xs0[i] + dh * tdir[i];
    xt1[i] = xs0[i] - dh * tdir[i];
    rs0[i] = 0.;
    hs0[i] = 1. - sdir[i] - tdir[i];
  }
  fe[0] = impl_func(xs0);
  fp0 = ivga.iat * impl_func(xs0);
  fs1 = ivga.iat * impl_func(xs1);
  fs2 = ivga.iat * impl_func(xs2);
  ft1 = ivga.iat * impl_func(xt1);
  ft2 = ivga.iat * impl_func(xt2);

  /* initial residue (-grad f) and diagonal hessian (in 2D) */
  dfs = -0.5 * (fs2 - fs1) / dh;
  dft = -0.5 * (ft2 - ft1) / dh;
  d2fs = (fs1 + fs2 - 2. * fp0) / (dh * dh);
  d2ft = (ft1 + ft2 - 2. * fp0) / (dh * dh);
  if (d2fs <= 0. || d2ft <= 0.) /* hes must be positive-definite */
    d2fs = d2ft = 1.;

  mcd = del0 = 0.;
  for (i = 0; i < NDIM; i++) {
    res[i] = rs0[i] + dfs * ivga.igs * sdir[i] + dft * ivga.igt * tdir[i];
    hes[i] = hs0[i] + d2fs * sdir[i] + d2ft * tdir[i];
    pcrs[i] = res[i] / hes[i]; /* residue/preconditioner */
    mcd += pcrs[i] * pcrs[i];
    del0 += res[i] * pcrs[i]; /* initial residue */
  }

  /* distance ss0 to the cell boundary */
  mcd = sqrt(mcd + EPS_NOT0);
  for (i = 0; i < NDIM; i++) {
    cndr[i] = pcrs[i]; /* conjugated direction */
    nmdr[i] = cndr[i] / mcd; /* unit conj. direction */
    d1 = SGN0P(nmdr[i]);
    d2 = fabs(nmdr[i]) + EPS_NOT0;
    if (d2 < EPS_R)
      ss[i] = 1000. * h0;
    else {
      a1 = (x0[i] - xs0[i]) / (d1 * d2);
      a2 = (x0[i] + h0 - xs0[i]) / (d1 * d2);
      ss[i] = MAX(a1, a2);
    }
  }
  ss0 = MIN(ss[0], ss[1]);
  ss0 = MIN(ss0, ss[2]);
  for (i = 0; i < NDIM; i++)
    xs1[i] = xs0[i] + ss0 * nmdr[i];
  fe[1] = impl_func(xs1);

  delnew = del0;
  not_conv = 1;
  iter = k = ipt = 0;
  while (not_conv && iter < max_iter) { /* iterative loop */
    /* DEBUG 1 */

    xfsa = vofi_get_segment_min(
        impl_func, fe, xs0, nmdr, ss0, ivga.iat, max_iter_line);

    for (i = 0; i < NDIM; i++)
      xs0[i] = xfsa.xval[i];
    fe[0] = xfsa.fval;
    fp0 = ivga.iat * xfsa.fval;
    if (xfsa.iat)
      not_conv = 0; /* got a sign change ! */
    else {
      for (i = 0; i < NDIM; i++) { /* discrete (- grad f) */
        xs2[i] = xs0[i] + dh * sdir[i];
        xs1[i] = xs0[i] - dh * sdir[i];
        xt2[i] = xs0[i] + dh * tdir[i];
        xt1[i] = xs0[i] - dh * tdir[i];
      }
      ss0 = xfsa.sval;
      /* DEBUG 2 */

      fs1 = ivga.iat * impl_func(xs1);
      fs2 = ivga.iat * impl_func(xs2);
      ft1 = ivga.iat * impl_func(xt1);
      ft2 = ivga.iat * impl_func(xt2);
      dfs = -0.5 * (fs2 - fs1) / dh;
      dft = -0.5 * (ft2 - ft1) / dh;
      d2fs = (fs1 + fs2 - 2. * fp0) / (dh * dh);
      d2ft = (ft1 + ft2 - 2. * fp0) / (dh * dh);
      if (d2fs <= 0. || d2ft <= 0.) d2fs = d2ft = 1.;
      delold = delnew;
      delmid = delnew = 0.;
      for (i = 0; i <= 2; i++) {
        res[i] = rs0[i] + dfs * sdir[i] + dft * tdir[i];
        delmid += res[i] * pcrs[i];
        hes[i] = hs0[i] + d2fs * sdir[i] + d2ft * tdir[i];
        pcrs[i] = res[i] / hes[i];
        delnew += res[i] * pcrs[i];
      }

      beta = (delnew - delmid) / delold;
      k++;
      /* DEBUG 3 */

      if (k == 2 || beta <= 0.) {
        beta = 0.;
        k = 0;
      }
      mcd = 0.;
      for (i = 0; i < NDIM; i++) {
        cndr[i] = pcrs[i] + beta * cndr[i];
        mcd += cndr[i] * cndr[i];
      }
      mcd = sqrt(mcd + EPS_NOT0);
      for (i = 0; i < NDIM; i++) {
        nmdr[i] = cndr[i] / mcd; /* unit conjugate direction */
        d1 = SGN0P(nmdr[i]);
        d2 = fabs(nmdr[i]) + EPS_NOT0;
        a1 = (x0[i] - xs0[i]) / (d1 * d2);
        a2 = (x0[i] + h0 - xs0[i]) / (d1 * d2);
        ss[i] = MAX(a1, a2);
      }
      ss1 = MIN(ss[0], ss[1]);
      ss1 = MIN(ss1, ss[2]);
      ss0 = MIN(1.2 * ss0, ss1);

      /* convergence criterion on residue or minimum on boundary */
      if (delnew < eps2 * del0 || ss0 < EPS_R) {
        not_conv = 0;
        /* DEBUG 4 */

      } else {
        for (i = 0; i < NDIM; i++)
          xs1[i] = xs0[i] + ss0 * nmdr[i];
        fe[1] = impl_func(xs1);
        iss = 0; /* bracket the minimum */
        while (ivga.iat * fe[1] < fp0 && iss < 3 && ss0 < ss1) {
          ss0 = MIN(3. * ss0, ss1);
          if (iss == 2) ss0 = ss1;
          for (i = 0; i < NDIM; i++)
            xs1[i] = xs0[i] + ss0 * nmdr[i];
          fe[1] = impl_func(xs1);
          iss++;
        }
      }
      iter++;
    }
  }

  return xfsa;
}
