#include <cmath>

#include "vect.hpp"

using Scal = double;
using Vect = geom::GVect<Scal, 3>;

#define MAX(a, b) (std::max(a, b))
#define MIN(a, b) (std::min(a, b))


void Clip(Scal& a, Scal l, Scal u) {
  a = std::max(l, std::min(u, a));
}

void Clip(Scal& a) {
  Clip(a, 0., 1.);
};


/**
 * gfs_line_area:
 * @m: normal to the line.
 * @alpha: line constant.
 *
 * Returns: the area of the fraction of a cell lying under the line
 * (@m,@alpha).
 */
Scal gfs_line_area (const Vect& m, Scal alpha)
{
  Vect n;
  Scal alpha1, a, v, area;

  //g_return_val_if_fail (m != NULL, 0.);

  n = m;
  alpha1 = alpha;
  if (n[0] < 0.) {
    alpha1 -= n[0];
    n[0] = - n[0];
  }
  if (n[1] < 0.) {
    alpha1 -= n[1];
    n[1] = - n[1];
  }

  if (alpha1 <= 0.)
    return 0.;

  if (alpha1 >= n[0] + n[1])
    return 1.;

  if (n[0] == 0.)
    area = alpha1/n[1];
  else if (n[1] == 0.)
    area = alpha1/n[0];
  else {
    v = alpha1*alpha1;

    a = alpha1 - n[0];
    if (a > 0.)
      v -= a*a;
    
    a = alpha1 - n[1];
    if (a > 0.)
      v -= a*a;

    area = v/(2.*n[0]*n[1]);
  }

  Clip(area);

  return area;
}

/**
 * gfs_line_alpha:
 * @m: a #Vect.
 * @c: a volume fraction.
 *
 * Returns: the value @alpha such that the area of a square cell
 * lying under the line defined by @m.@x = @alpha is equal to @c. 
 */
Scal gfs_line_alpha (const Vect& m, Scal c)
{
  Scal alpha, m1, m2, v1;

  //g_return_val_if_fail (m != NULL, 0.);
  //g_return_val_if_fail (c >= 0. && c <= 1., 0.);
  
  m1 = fabs (m[0]); m2 = fabs (m[1]);
  if (m1 > m2) {
    v1 = m1; m1 = m2; m2 = v1;
  }
  
  v1 = m1/2.;
  if (c <= v1/m2)
    alpha = sqrt (2.*c*m1*m2);
  else if (c <= 1. - v1/m2)
    alpha = c*m2 + v1;
  else
    alpha = m1 + m2 - sqrt (2.*m1*m2*(1. - c));

  if (m[0] < 0.)
    alpha += m[0];
  if (m[1] < 0.)
    alpha += m[1];

  return alpha;
}

#define EPS 1e-4

/**
 * gfs_line_center:
 * @m: normal to the line.
 * @alpha: line constant.
 * @a: area of cell fraction.
 * @p: a #Vect.
 *
 * Fills @p with the position of the center of mass of the fraction of
 * a square cell lying under the line (@m,@alpha).
 */
void gfs_line_center (const Vect& m, Scal alpha, Scal a, Vect& p)
{
  Vect n;
  Scal b;

  //g_return_if_fail (m != NULL);
  //g_return_if_fail (p != NULL);

  n = m;
  if (n[0] < 0.) {
    alpha -= n[0];
    n[0] = - n[0];
  }
  if (n[1] < 0.) {
    alpha -= n[1];
    n[1] = - n[1];
  }

  p[2] = 0.;
  if (alpha <= 0.) {
    p[0] = p[1] = 0.;
    return;
  }

  if (alpha >= n[0] + n[1]) {
    p[0] = p[1] = 0.5;
    return;
  }

  //g_return_if_fail (a > 0. && a < 1.);

  if (n[0] < EPS) {
    p[0] = 0.5;
    p[1] = m[1] < 0. ? 1. - a/2. : a/2.;
    return;
  }

  if (n[1] < EPS) {
    p[1] = 0.5;
    p[0] = m[0] < 0. ? 1. - a/2. : a/2.;
    return;
  }

  p[0] = p[1] = alpha*alpha*alpha;

  b = alpha - n[0];
  if (b > 0.) {
    p[0] -= b*b*(alpha + 2.*n[0]);
    p[1] -= b*b*b;
  }

  b = alpha - n[1];
  if (b > 0.) {
    p[1] -= b*b*(alpha + 2.*n[1]);
    p[0] -= b*b*b;
  }
  
  p[0] /= 6.*n[0]*n[0]*n[1]*a;
  p[1] /= 6.*n[0]*n[1]*n[1]*a;

  if (m[0] < 0.)
    p[0] = 1. - p[0];
  if (m[1] < 0.)
    p[1] = 1. - p[1];
}

/**
 * gfs_line_area_center:
 * @m: normal to the line.
 * @alpha: line constant.
 * @p: a #Vect.
 *
 * Fills @p with the position of the center of area of the fraction of
 * a square cell lying under the line (@m,@alpha).
 *
 * Returns: the length of the facet.
 */
Scal gfs_line_area_center (const Vect& m, Scal alpha, Vect& p)
{
  Vect n;

  //g_return_val_if_fail (m != NULL, 0.);
  //g_return_val_if_fail (p != NULL, 0.);

  n = m;
  if (n[0] < 0.) {
    alpha -= n[0];
    n[0] = - n[0];
  }
  if (n[1] < 0.) {
    alpha -= n[1];
    n[1] = - n[1];
  }

  p[2] = 0.;
  if (alpha <= 0. || alpha >= n[0] + n[1]) {
    p[0] = p[1] = 0.;
    return 0.;
  }

  if (n[0] < EPS) {
    p[0] = 0.5;
    p[1] = m[1] < 0. ? 1. - alpha : alpha;
    return 1.;
  }

  if (n[1] < EPS) {
    p[1] = 0.5;
    p[0] = m[0] < 0. ? 1. - alpha : alpha;
    return 1.;
  }

  p[0] = p[1] = 0.;

  if (alpha >= n[0]) {
    p[0] += 1.;
    p[1] += (alpha - n[0])/n[1];
  }
  else
    p[0] += alpha/n[0];

  Scal ax = p[0], ay = p[1];
  if (alpha >= n[1]) {
    p[1] += 1.;
    ay -= 1.;
    p[0] += (alpha - n[1])/n[0];
    ax -= (alpha - n[1])/n[0];
  }
  else {
    p[1] += alpha/n[1];
    ay -= alpha/n[1];
  }

  p[0] /= 2.;
  p[1] /= 2.;

  Clip(p[0]);
  Clip(p[1]);

  if (m[0] < 0.)
    p[0] = 1. - p[0];
  if (m[1] < 0.)
    p[1] = 1. - p[1];

  return sqrt (ax*ax + ay*ay);
}

/**
 * gfs_plane_volume:
 * @m: normal to the plane.
 * @alpha: plane constant.
 *
 * Returns: the volume of a cell lying under the plane (@m,@alpha).
 */
Scal gfs_plane_volume (const Vect& m, Scal alpha)
{
  //g_return_val_if_fail (m != NULL, 0.);

  Scal al = alpha + MAX(0., -m[0]) + MAX(0., -m[1]) + MAX(0., -m[2]);
  if (al <= 0.)
    return 0.;
  Scal tmp = fabs(m[0]) + fabs(m[1]) + fabs(m[2]);
  if (al >= tmp)
    return 1.;
  //g_assert (tmp > 0.);
  Scal n1 = fabs(m[0])/tmp;
  Scal n2 = fabs(m[1])/tmp;
  Scal n3 = fabs(m[2])/tmp;
  al = MAX(0., MIN(1., al/tmp));
  Scal al0 = MIN(al, 1. - al);
  Scal b1 = MIN(n1*1, n2);
  Scal b3 = MAX(n1*1, n2);
  Scal b2 = n3;
  if (b2 < b1) {
    tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  else if (b2 > b3) {
    tmp = b3;
    b3 = b2;
    b2 = tmp;
  }
  Scal b12 = b1 + b2;
  Scal bm = MIN(b12, b3);
  Scal pr = MAX(6.*b1*b2*b3, 1e-50);
  if (al0 < b1)
    tmp = al0*al0*al0/pr;
  else if (al0 < b2)
    tmp = 0.5*al0*(al0 - b1)/(b2*b3) +  b1*b1*b1/pr;
  else if (al0 < bm)
    tmp = (al0*al0*(3.*b12 - al0) + b1*b1*(b1 - 3.*al0) + b2*b2*(b2 - 3.*al0))/pr;
  else if (b12 < b3)
    tmp = (al0 - 0.5*bm)/b3;
  else
    tmp = (al0*al0*(3. - 2.*al0) + b1*b1*(b1 - 3.*al0) + 
	   b2*b2*(b2 - 3.*al0) + b3*b3*(b3 - 3.*al0))/pr;

  Scal volume = al <= 0.5 ? tmp : 1. - tmp;
  Clip(volume);
  return volume;
}

/**
 * gfs_plane_alpha:
 * @m: a #Vect.
 * @c: a volume fraction.
 *
 * Returns: the value @alpha such that the volume of a cubic cell
 * lying under the plane defined by @m.@x = @alpha is equal to @c. 
 */
Scal gfs_plane_alpha (const Vect& m, Scal c)
{
  Scal alpha;
  Vect n;

  //g_return_val_if_fail (m != NULL, 0.);
  //g_return_val_if_fail (c >= 0. && c <= 1., 0.);

  n[0] = fabs (m[0]); n[1] = fabs (m[1]); n[2] = fabs (m[2]);

  Scal m1, m2, m3;
  m1 = MIN(n[0], n[1]);
  m3 = MAX(n[0], n[1]);
  m2 = n[2];
  if (m2 < m1) {
    Scal tmp = m1;
    m1 = m2;
    m2 = tmp;
  }
  else if (m2 > m3) {
    Scal tmp = m3;
    m3 = m2;
    m2 = tmp;
  }
  Scal m12 = m1 + m2;
  Scal pr = MAX(6.*m1*m2*m3, 1e-50);
  Scal V1 = m1*m1*m1/pr;
  Scal V2 = V1 + (m2 - m1)/(2.*m3), V3;
  Scal mm;
  if (m3 < m12) {
    mm = m3;
    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;
  }
  else {
    mm = m12;
    V3 = mm/(2.*m3);
  }

  Scal ch = MIN(c, 1. - c);
  if (ch < V1)
    alpha = pow (pr*ch, 1./3.);
  else if (ch < V2)
    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;
  else if (ch < V3) {
    Scal p = 2.*m1*m2;
    Scal q = 3.*m1*m2*(m12 - 2.*m3*ch)/2.;
    Scal p12 = sqrt (p);
    Scal teta = acos(q/(p*p12))/3.;
    Scal cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;
  }
  else if (m12 < m3)
    alpha = m3*ch + mm/2.;
  else {
    Scal p = m1*(m2 + m3) + m2*m3 - 1./4.;
    Scal q = 3.*m1*m2*m3*(1./2. - ch)/2.;
    Scal p12 = sqrt(p);
    Scal teta = acos(q/(p*p12))/3.;
    Scal cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }
  if (c > 1./2.) alpha = 1. - alpha;

  if (m[0] < 0.)
    alpha += m[0];
  if (m[1] < 0.)
    alpha += m[1];
  if (m[2] < 0.)
    alpha += m[2];

  return alpha;
}

/**
 * gfs_plane_center:
 * @m: normal to the plane.
 * @alpha: plane constant.
 * @a: volume of cell fraction.
 * @p: a #Vect.
 *
 * Fills @p with the position of the center of mass of the fraction of
 * a cubic cell lying under the plane (@m,@alpha).
 */
void gfs_plane_center (const Vect& m, Scal alpha, Scal a, Vect& p)
{
  Vect n;
  Scal b, amax;

  //g_return_if_fail (m != NULL);
  //g_return_if_fail (p != NULL);
  //g_return_if_fail (a >= 0. && a <= 1.);

  if (fabs (m[0]) < EPS) {
    Vect q;
    n[0] = m[1];
    n[1] = m[2];
    gfs_line_center (n, alpha, a, q);
    p[0] = 0.5;
    p[1] = q[0];
    p[2] = q[1];
    return;
  }
  if (fabs (m[1]) < EPS) {
    Vect q;
    n[0] = m[2];
    n[1] = m[0];
    gfs_line_center (n, alpha, a, q);
    p[0] = q[1];
    p[1] = 0.5;
    p[2] = q[0];
    return;
  }
  if (fabs (m[2]) < EPS) {
    gfs_line_center (m, alpha, a, p);
    p[2] = 0.5;
    return;
  }

  n = m;
  if (n[0] < 0.) {
    alpha -= n[0];
    n[0] = - n[0];
  }
  if (n[1] < 0.) {
    alpha -= n[1];
    n[1] = - n[1];
  }  
  if (n[2] < 0.) {
    alpha -= n[2];
    n[2] = - n[2];
  }  

  if (alpha <= 0. || a == 0.) {
    p[0] = p[1] = p[2] = 0.;
    return;
  }

  if (alpha >= n[0] + n[1] + n[2] || a == 1.) {
    p[0] = p[1] = p[2] = 0.5;
    return;
  }

  amax = n[0] + n[1] + n[2];
  p[0] = p[1] = p[2] = alpha*alpha*alpha*alpha;

  b = alpha - n[0];
  if (b > 0.) {
    p[0] -= b*b*b*(3.*n[0] + alpha);
    p[1] -= b*b*b*b;
    p[2] -= b*b*b*b;
  }
  b = alpha - n[1];
  if (b > 0.) {
    p[1] -= b*b*b*(3.*n[1] + alpha);
    p[0] -= b*b*b*b;
    p[2] -= b*b*b*b;
  }
  b = alpha - n[2];
  if (b > 0.) {
    p[2] -= b*b*b*(3.*n[2] + alpha);
    p[0] -= b*b*b*b;
    p[1] -= b*b*b*b;
  }

  amax = alpha - amax;
  b = amax + n[0];
  if (b > 0.) {
    p[1] += b*b*b*(3.*n[1] + alpha - n[2]);
    p[2] += b*b*b*(3.*n[2] + alpha - n[1]);
    p[0] += b*b*b*b;
  }
  b = amax + n[1];
  if (b > 0.) {
    p[0] += b*b*b*(3.*n[0] + alpha - n[2]);
    p[2] += b*b*b*(3.*n[2] + alpha - n[0]);
    p[1] += b*b*b*b;
  }
  b = amax + n[2];
  if (b > 0.) {
    p[0] += b*b*b*(3.*n[0] + alpha - n[1]);
    p[1] += b*b*b*(3.*n[1] + alpha - n[0]);
    p[2] += b*b*b*b;
  }

  b = 24.*n[0]*n[1]*n[2]*a;
  p[0] /= b*n[0]; p[1] /= b*n[1]; p[2] /= b*n[2];

  if (m[0] < 0.) p[0] = 1. - p[0];
  if (m[1] < 0.) p[1] = 1. - p[1];
  if (m[2] < 0.) p[2] = 1. - p[2];
}

/**
 * gfs_plane_area_center:
 * @m: normal to the plane.
 * @alpha: plane constant.
 * @p: a #Vect.
 *
 * Fills @p with the position of the center of area of the fraction of
 * a cubic cell lying under the plane (@m,@alpha).
 *
 * Returns: the area of the facet.
 */
Scal gfs_plane_area_center (const Vect& m, Scal alpha, Vect& p)
{
  //g_return_val_if_fail (m != NULL, 0.);
  //g_return_val_if_fail (p != NULL, 0.);

  if (fabs (m[0]) < EPS) {
    Vect n, q;
    n[0] = m[1];
    n[1] = m[2];
    Scal area = gfs_line_area_center (n, alpha, q);
    p[0] = 0.5;
    p[1] = q[0];
    p[2] = q[1];
    return area;
  }
  if (fabs (m[1]) < EPS) {
    Vect n, q;
    n[0] = m[2];
    n[1] = m[0];
    Scal area = gfs_line_area_center (n, alpha, q);
    p[0] = q[1];
    p[1] = 0.5;
    p[2] = q[0];
    return area;
  }
  if (fabs (m[2]) < EPS) {
    Scal area = gfs_line_area_center (m, alpha, p);
    p[2] = 0.5;
    return area;
  }

  Vect n = m;
  if (n[0] < 0.) {
    alpha -= n[0];
    n[0] = - n[0];
  }
  if (n[1] < 0.) {
    alpha -= n[1];
    n[1] = - n[1];
  }  
  if (n[2] < 0.) {
    alpha -= n[2];
    n[2] = - n[2];
  }

  Scal amax = n[0] + n[1] + n[2];
  if (alpha <= 0. || alpha >= amax) {
    p[0] = p[1] = p[2] = 0.;
    return 0.;
  }

  Scal area = alpha*alpha;
  p[0] = p[1] = p[2] = area*alpha;

  Scal b = alpha - n[0];
  if (b > 0.) {
    area -= b*b;
    p[0] -= b*b*(2.*n[0] + alpha);
    p[1] -= b*b*b;
    p[2] -= b*b*b;
  }
  b = alpha - n[1];
  if (b > 0.) {
    area -= b*b;
    p[1] -= b*b*(2.*n[1] + alpha);
    p[0] -= b*b*b;
    p[2] -= b*b*b;
  }
  b = alpha - n[2];
  if (b > 0.) {
    area -= b*b;
    p[2] -= b*b*(2.*n[2] + alpha);
    p[0] -= b*b*b;
    p[1] -= b*b*b;
  }

  amax = alpha - amax;
  b = amax + n[0];
  if (b > 0.) {
    area += b*b;
    p[1] += b*b*(2.*n[1] + alpha - n[2]);
    p[2] += b*b*(2.*n[2] + alpha - n[1]);
    p[0] += b*b*b;
  }
  b = amax + n[1];
  if (b > 0.) {
    area += b*b;
    p[0] += b*b*(2.*n[0] + alpha - n[2]);
    p[2] += b*b*(2.*n[2] + alpha - n[0]);
    p[1] += b*b*b;
  }
  b = amax + n[2];
  if (b > 0.) {
    area += b*b;
    p[0] += b*b*(2.*n[0] + alpha - n[1]);
    p[1] += b*b*(2.*n[1] + alpha - n[0]);
    p[2] += b*b*b;
  }

  area *= 3.;
  p[0] /= area*n[0];
  p[1] /= area*n[1];
  p[2] /= area*n[2];

  Clip(p[0]);
  Clip(p[1]);
  Clip(p[2]);

  if (m[0] < 0.) p[0] = 1. - p[0];
  if (m[1] < 0.) p[1] = 1. - p[1];
  if (m[2] < 0.) p[2] = 1. - p[2];

  return area*sqrt (1./(n[0]*n[0]*n[1]*n[1]) + 1./(n[0]*n[0]*n[2]*n[2]) + 1./(n[2]*n[2]*n[1]*n[1]))/6.;
}

