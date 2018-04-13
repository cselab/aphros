#include <cmath>

#include "reconst.h"

#define MAX(a, b) (std::max(a, b))
#define MIN(a, b) (std::min(a, b))
#define EPS 1e-4


void Clip(Scal& a, Scal l, Scal u) {
  a = std::max(l, std::min(u, a));
}

void Clip(Scal& a) {
  Clip(a, 0., 1.);
};


Scal gfs_line_s (const Vect& m, Scal a)
{
  Vect n = m;
  Scal a1 = a;
  if (n[0] < 0.) {
    a1 -= n[0];
    n[0] = - n[0];
  }
  if (n[1] < 0.) {
    a1 -= n[1];
    n[1] = - n[1];
  }

  if (a1 <= 0.)
    return 0.;

  if (a1 >= n[0] + n[1])
    return 1.;

  Scal s;
  if (n[0] == 0.)
    s = a1 / n[1];
  else if (n[1] == 0.)
    s = a1 / n[0];
  else {
    Scal v = a1 * a1;

    Scal d = a1 - n[0];
    if (a > 0.)
      v -= d * d;
    
    a = a1 - n[1];
    if (a > 0.)
      v -= d * d;

    s = v / (2. * n[0] * n[1]);
  }

  Clip(s);

  return s;
}

Scal gfs_line_a (const Vect& m, Scal c)
{
  Scal a, m1, m2, v1;

  //g_return_val_if_fail (m != NULL, 0.);
  //g_return_val_if_fail (c >= 0. && c <= 1., 0.);
  
  m1 = fabs (m[0]); m2 = fabs (m[1]);
  if (m1 > m2) {
    v1 = m1; m1 = m2; m2 = v1;
  }
  
  v1 = m1/2.;
  if (c <= v1/m2)
    a = sqrt (2.*c*m1*m2);
  else if (c <= 1. - v1/m2)
    a = c*m2 + v1;
  else
    a = m1 + m2 - sqrt (2.*m1*m2*(1. - c));

  if (m[0] < 0.)
    a += m[0];
  if (m[1] < 0.)
    a += m[1];

  return a;
}

#define EPS 1e-4

void gfs_line_c (const Vect& m, Scal a, Scal s, Vect& p)
{
  Vect n;
  Scal b;

  //g_return_if_fail (m != NULL);
  //g_return_if_fail (p != NULL);

  n = m;
  if (n[0] < 0.) {
    a -= n[0];
    n[0] = - n[0];
  }
  if (n[1] < 0.) {
    a -= n[1];
    n[1] = - n[1];
  }

  p[2] = 0.;
  if (a <= 0.) {
    p[0] = p[1] = 0.;
    return;
  }

  if (a >= n[0] + n[1]) {
    p[0] = p[1] = 0.5;
    return;
  }

  //g_return_if_fail (s > 0. && s < 1.);

  if (n[0] < EPS) {
    p[0] = 0.5;
    p[1] = m[1] < 0. ? 1. - s/2. : s/2.;
    return;
  }

  if (n[1] < EPS) {
    p[1] = 0.5;
    p[0] = m[0] < 0. ? 1. - s/2. : s/2.;
    return;
  }

  p[0] = p[1] = a*a*a;

  b = a - n[0];
  if (b > 0.) {
    p[0] -= b*b*(a + 2.*n[0]);
    p[1] -= b*b*b;
  }

  b = a - n[1];
  if (b > 0.) {
    p[1] -= b*b*(a + 2.*n[1]);
    p[0] -= b*b*b;
  }
  
  p[0] /= 6.*n[0]*n[0]*n[1]*s;
  p[1] /= 6.*n[0]*n[1]*n[1]*s;

  if (m[0] < 0.)
    p[0] = 1. - p[0];
  if (m[1] < 0.)
    p[1] = 1. - p[1];
}

Scal gfs_line_sc (const Vect& m, Scal a, Vect& p)
{
  Vect n;

  //g_return_val_if_fail (m != NULL, 0.);
  //g_return_val_if_fail (p != NULL, 0.);

  n = m;
  if (n[0] < 0.) {
    a -= n[0];
    n[0] = - n[0];
  }
  if (n[1] < 0.) {
    a -= n[1];
    n[1] = - n[1];
  }

  p[2] = 0.;
  if (a <= 0. || a >= n[0] + n[1]) {
    p[0] = p[1] = 0.;
    return 0.;
  }

  if (n[0] < EPS) {
    p[0] = 0.5;
    p[1] = m[1] < 0. ? 1. - a : a;
    return 1.;
  }

  if (n[1] < EPS) {
    p[1] = 0.5;
    p[0] = m[0] < 0. ? 1. - a : a;
    return 1.;
  }

  p[0] = p[1] = 0.;

  if (a >= n[0]) {
    p[0] += 1.;
    p[1] += (a - n[0])/n[1];
  }
  else
    p[0] += a/n[0];

  Scal ax = p[0], ay = p[1];
  if (a >= n[1]) {
    p[1] += 1.;
    ay -= 1.;
    p[0] += (a - n[1])/n[0];
    ax -= (a - n[1])/n[0];
  }
  else {
    p[1] += a/n[1];
    ay -= a/n[1];
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

Scal gfs_plane_volume (const Vect& m, Scal a)
{
  //g_return_val_if_fail (m != NULL, 0.);

  Scal al = a + MAX(0., -m[0]) + MAX(0., -m[1]) + MAX(0., -m[2]);
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

Scal gfs_plane_a (const Vect& m, Scal c)
{
  Scal a;
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
    a = pow (pr*ch, 1./3.);
  else if (ch < V2)
    a = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;
  else if (ch < V3) {
    Scal p = 2.*m1*m2;
    Scal q = 3.*m1*m2*(m12 - 2.*m3*ch)/2.;
    Scal p12 = sqrt (p);
    Scal teta = acos(q/(p*p12))/3.;
    Scal cs = cos(teta);
    a = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;
  }
  else if (m12 < m3)
    a = m3*ch + mm/2.;
  else {
    Scal p = m1*(m2 + m3) + m2*m3 - 1./4.;
    Scal q = 3.*m1*m2*m3*(1./2. - ch)/2.;
    Scal p12 = sqrt(p);
    Scal teta = acos(q/(p*p12))/3.;
    Scal cs = cos(teta);
    a = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }
  if (c > 1./2.) a = 1. - a;

  if (m[0] < 0.)
    a += m[0];
  if (m[1] < 0.)
    a += m[1];
  if (m[2] < 0.)
    a += m[2];

  return a;
}

void gfs_plane_c (const Vect& m, Scal a, Scal s, Vect& p)
{
  Vect n;
  Scal b, amax;

  //g_return_if_fail (m != NULL);
  //g_return_if_fail (p != NULL);
  //g_return_if_fail (s >= 0. && s <= 1.);

  if (fabs (m[0]) < EPS) {
    Vect q;
    n[0] = m[1];
    n[1] = m[2];
    gfs_line_c (n, a, s, q);
    p[0] = 0.5;
    p[1] = q[0];
    p[2] = q[1];
    return;
  }
  if (fabs (m[1]) < EPS) {
    Vect q;
    n[0] = m[2];
    n[1] = m[0];
    gfs_line_c (n, a, s, q);
    p[0] = q[1];
    p[1] = 0.5;
    p[2] = q[0];
    return;
  }
  if (fabs (m[2]) < EPS) {
    gfs_line_c (m, a, s, p);
    p[2] = 0.5;
    return;
  }

  n = m;
  if (n[0] < 0.) {
    a -= n[0];
    n[0] = - n[0];
  }
  if (n[1] < 0.) {
    a -= n[1];
    n[1] = - n[1];
  }  
  if (n[2] < 0.) {
    a -= n[2];
    n[2] = - n[2];
  }  

  if (a <= 0. || s == 0.) {
    p[0] = p[1] = p[2] = 0.;
    return;
  }

  if (a >= n[0] + n[1] + n[2] || s == 1.) {
    p[0] = p[1] = p[2] = 0.5;
    return;
  }

  amax = n[0] + n[1] + n[2];
  p[0] = p[1] = p[2] = a*a*a*a;

  b = a - n[0];
  if (b > 0.) {
    p[0] -= b*b*b*(3.*n[0] + a);
    p[1] -= b*b*b*b;
    p[2] -= b*b*b*b;
  }
  b = a - n[1];
  if (b > 0.) {
    p[1] -= b*b*b*(3.*n[1] + a);
    p[0] -= b*b*b*b;
    p[2] -= b*b*b*b;
  }
  b = a - n[2];
  if (b > 0.) {
    p[2] -= b*b*b*(3.*n[2] + a);
    p[0] -= b*b*b*b;
    p[1] -= b*b*b*b;
  }

  amax = a - amax;
  b = amax + n[0];
  if (b > 0.) {
    p[1] += b*b*b*(3.*n[1] + a - n[2]);
    p[2] += b*b*b*(3.*n[2] + a - n[1]);
    p[0] += b*b*b*b;
  }
  b = amax + n[1];
  if (b > 0.) {
    p[0] += b*b*b*(3.*n[0] + a - n[2]);
    p[2] += b*b*b*(3.*n[2] + a - n[0]);
    p[1] += b*b*b*b;
  }
  b = amax + n[2];
  if (b > 0.) {
    p[0] += b*b*b*(3.*n[0] + a - n[1]);
    p[1] += b*b*b*(3.*n[1] + a - n[0]);
    p[2] += b*b*b*b;
  }

  b = 24.*n[0]*n[1]*n[2]*s;
  p[0] /= b*n[0]; p[1] /= b*n[1]; p[2] /= b*n[2];

  if (m[0] < 0.) p[0] = 1. - p[0];
  if (m[1] < 0.) p[1] = 1. - p[1];
  if (m[2] < 0.) p[2] = 1. - p[2];
}

Scal gfs_plane_sc (const Vect& m, Scal a, Vect& p)
{
  //g_return_val_if_fail (m != NULL, 0.);
  //g_return_val_if_fail (p != NULL, 0.);

  if (fabs (m[0]) < EPS) {
    Vect n, q;
    n[0] = m[1];
    n[1] = m[2];
    Scal s = gfs_line_sc (n, a, q);
    p[0] = 0.5;
    p[1] = q[0];
    p[2] = q[1];
    return s;
  }
  if (fabs (m[1]) < EPS) {
    Vect n, q;
    n[0] = m[2];
    n[1] = m[0];
    Scal s = gfs_line_sc (n, a, q);
    p[0] = q[1];
    p[1] = 0.5;
    p[2] = q[0];
    return s;
  }
  if (fabs (m[2]) < EPS) {
    Scal s = gfs_line_sc (m, a, p);
    p[2] = 0.5;
    return s;
  }

  Vect n = m;
  if (n[0] < 0.) {
    a -= n[0];
    n[0] = - n[0];
  }
  if (n[1] < 0.) {
    a -= n[1];
    n[1] = - n[1];
  }  
  if (n[2] < 0.) {
    a -= n[2];
    n[2] = - n[2];
  }

  Scal amax = n[0] + n[1] + n[2];
  if (a <= 0. || a >= amax) {
    p[0] = p[1] = p[2] = 0.;
    return 0.;
  }

  Scal s = a*a;
  p[0] = p[1] = p[2] = s*a;

  Scal b = a - n[0];
  if (b > 0.) {
    s -= b*b;
    p[0] -= b*b*(2.*n[0] + a);
    p[1] -= b*b*b;
    p[2] -= b*b*b;
  }
  b = a - n[1];
  if (b > 0.) {
    s -= b*b;
    p[1] -= b*b*(2.*n[1] + a);
    p[0] -= b*b*b;
    p[2] -= b*b*b;
  }
  b = a - n[2];
  if (b > 0.) {
    s -= b*b;
    p[2] -= b*b*(2.*n[2] + a);
    p[0] -= b*b*b;
    p[1] -= b*b*b;
  }

  amax = a - amax;
  b = amax + n[0];
  if (b > 0.) {
    s += b*b;
    p[1] += b*b*(2.*n[1] + a - n[2]);
    p[2] += b*b*(2.*n[2] + a - n[1]);
    p[0] += b*b*b;
  }
  b = amax + n[1];
  if (b > 0.) {
    s += b*b;
    p[0] += b*b*(2.*n[0] + a - n[2]);
    p[2] += b*b*(2.*n[2] + a - n[0]);
    p[1] += b*b*b;
  }
  b = amax + n[2];
  if (b > 0.) {
    s += b*b;
    p[0] += b*b*(2.*n[0] + a - n[1]);
    p[1] += b*b*(2.*n[1] + a - n[0]);
    p[2] += b*b*b;
  }

  s *= 3.;
  p[0] /= s*n[0];
  p[1] /= s*n[1];
  p[2] /= s*n[2];

  Clip(p[0]);
  Clip(p[1]);
  Clip(p[2]);

  if (m[0] < 0.) p[0] = 1. - p[0];
  if (m[1] < 0.) p[1] = 1. - p[1];
  if (m[2] < 0.) p[2] = 1. - p[2];

  return s*sqrt (1./(n[0]*n[0]*n[1]*n[1]) + 1./(n[0]*n[0]*n[2]*n[2]) + 1./(n[2]*n[2]*n[1]*n[1]))/6.;
}

