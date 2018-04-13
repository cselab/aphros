/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2011 National Institute of Water and Atmospheric Research
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.  
 */
/*! \file
 * \brief Volume-Of-Fluid tracers.
 */

#include <math.h>
#include <stdlib.h>
#include "vof.h"
#include "variable.h"
#include "adaptive.h"
#include "graphic.h"

#define THRESHOLD(c) {if ((c) < 0.) c = 0.; else if ((c) > 1.) c = 1.;}

/**
 * gfs_line_area:
 * @m: normal to the line.
 * @alpha: line constant.
 *
 * Returns: the area of the fraction of a cell lying under the line
 * (@m,@alpha).
 */
gdouble gfs_line_area (const FttVector * m, gdouble alpha)
{
  FttVector n;
  gdouble alpha1, a, v, area;

  g_return_val_if_fail (m != NULL, 0.);

  n = *m;
  alpha1 = alpha;
  if (n.x < 0.) {
    alpha1 -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha1 -= n.y;
    n.y = - n.y;
  }

  if (alpha1 <= 0.)
    return 0.;

  if (alpha1 >= n.x + n.y)
    return 1.;

  if (n.x == 0.)
    area = alpha1/n.y;
  else if (n.y == 0.)
    area = alpha1/n.x;
  else {
    v = alpha1*alpha1;

    a = alpha1 - n.x;
    if (a > 0.)
      v -= a*a;
    
    a = alpha1 - n.y;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*n.x*n.y);
  }

  return CLAMP (area, 0., 1.);
}

/**
 * gfs_line_alpha:
 * @m: a #FttVector.
 * @c: a volume fraction.
 *
 * Returns: the value @alpha such that the area of a square cell
 * lying under the line defined by @m.@x = @alpha is equal to @c. 
 */
gdouble gfs_line_alpha (const FttVector * m, gdouble c)
{
  gdouble alpha, m1, m2, v1;

  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (c >= 0. && c <= 1., 0.);
  
  m1 = fabs (m->x); m2 = fabs (m->y);
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

  if (m->x < 0.)
    alpha += m->x;
  if (m->y < 0.)
    alpha += m->y;

  return alpha;
}

#define EPS 1e-4

/**
 * gfs_line_center:
 * @m: normal to the line.
 * @alpha: line constant.
 * @a: area of cell fraction.
 * @p: a #FttVector.
 *
 * Fills @p with the position of the center of mass of the fraction of
 * a square cell lying under the line (@m,@alpha).
 */
void gfs_line_center (const FttVector * m, gdouble alpha, gdouble a, FttVector * p)
{
  FttVector n;
  gdouble b;

  g_return_if_fail (m != NULL);
  g_return_if_fail (p != NULL);

  n = *m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }

  p->z = 0.;
  if (alpha <= 0.) {
    p->x = p->y = 0.;
    return;
  }

  if (alpha >= n.x + n.y) {
    p->x = p->y = 0.5;
    return;
  }

  g_return_if_fail (a > 0. && a < 1.);

  if (n.x < EPS) {
    p->x = 0.5;
    p->y = m->y < 0. ? 1. - a/2. : a/2.;
    return;
  }

  if (n.y < EPS) {
    p->y = 0.5;
    p->x = m->x < 0. ? 1. - a/2. : a/2.;
    return;
  }

  p->x = p->y = alpha*alpha*alpha;

  b = alpha - n.x;
  if (b > 0.) {
    p->x -= b*b*(alpha + 2.*n.x);
    p->y -= b*b*b;
  }

  b = alpha - n.y;
  if (b > 0.) {
    p->y -= b*b*(alpha + 2.*n.y);
    p->x -= b*b*b;
  }
  
  p->x /= 6.*n.x*n.x*n.y*a;
  p->y /= 6.*n.x*n.y*n.y*a;

  if (m->x < 0.)
    p->x = 1. - p->x;
  if (m->y < 0.)
    p->y = 1. - p->y;
}

/**
 * gfs_line_area_center:
 * @m: normal to the line.
 * @alpha: line constant.
 * @p: a #FttVector.
 *
 * Fills @p with the position of the center of area of the fraction of
 * a square cell lying under the line (@m,@alpha).
 *
 * Returns: the length of the facet.
 */
gdouble gfs_line_area_center (const FttVector * m, gdouble alpha, FttVector * p)
{
  FttVector n;

  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (p != NULL, 0.);

  n = *m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }

  p->z = 0.;
  if (alpha <= 0. || alpha >= n.x + n.y) {
    p->x = p->y = 0.;
    return 0.;
  }

  if (n.x < EPS) {
    p->x = 0.5;
    p->y = m->y < 0. ? 1. - alpha : alpha;
    return 1.;
  }

  if (n.y < EPS) {
    p->y = 0.5;
    p->x = m->x < 0. ? 1. - alpha : alpha;
    return 1.;
  }

  p->x = p->y = 0.;

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  gdouble ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

  p->x /= 2.;
  p->y /= 2.;

  THRESHOLD (p->x);
  THRESHOLD (p->y);

  if (m->x < 0.)
    p->x = 1. - p->x;
  if (m->y < 0.)
    p->y = 1. - p->y;

  return sqrt (ax*ax + ay*ay);
}

#if (!FTT_2D)
/**
 * gfs_plane_volume:
 * @m: normal to the plane.
 * @alpha: plane constant.
 *
 * Returns: the volume of a cell lying under the plane (@m,@alpha).
 */
gdouble gfs_plane_volume (const FttVector * m, gdouble alpha)
{
  g_return_val_if_fail (m != NULL, 0.);

  gdouble al = alpha + MAX(0., -m->x) + MAX(0., -m->y) + MAX(0., -m->z);
  if (al <= 0.)
    return 0.;
  gdouble tmp = fabs(m->x) + fabs(m->y) + fabs(m->z);
  if (al >= tmp)
    return 1.;
  g_assert (tmp > 0.);
  gdouble n1 = fabs(m->x)/tmp;
  gdouble n2 = fabs(m->y)/tmp;
  gdouble n3 = fabs(m->z)/tmp;
  al = MAX(0., MIN(1., al/tmp));
  gdouble al0 = MIN(al, 1. - al);
  gdouble b1 = MIN(n1*1, n2);
  gdouble b3 = MAX(n1*1, n2);
  gdouble b2 = n3;
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
  gdouble b12 = b1 + b2;
  gdouble bm = MIN(b12, b3);
  gdouble pr = MAX(6.*b1*b2*b3, 1e-50);
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

  gdouble volume = al <= 0.5 ? tmp : 1. - tmp;
  return CLAMP (volume, 0., 1.);
}

/**
 * gfs_plane_alpha:
 * @m: a #FttVector.
 * @c: a volume fraction.
 *
 * Returns: the value @alpha such that the volume of a cubic cell
 * lying under the plane defined by @m.@x = @alpha is equal to @c. 
 */
gdouble gfs_plane_alpha (const FttVector * m, gdouble c)
{
  gdouble alpha;
  FttVector n;

  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (c >= 0. && c <= 1., 0.);

  n.x = fabs (m->x); n.y = fabs (m->y); n.z = fabs (m->z);

  gdouble m1, m2, m3;
  m1 = MIN(n.x, n.y);
  m3 = MAX(n.x, n.y);
  m2 = n.z;
  if (m2 < m1) {
    gdouble tmp = m1;
    m1 = m2;
    m2 = tmp;
  }
  else if (m2 > m3) {
    gdouble tmp = m3;
    m3 = m2;
    m2 = tmp;
  }
  gdouble m12 = m1 + m2;
  gdouble pr = MAX(6.*m1*m2*m3, 1e-50);
  gdouble V1 = m1*m1*m1/pr;
  gdouble V2 = V1 + (m2 - m1)/(2.*m3), V3;
  gdouble mm;
  if (m3 < m12) {
    mm = m3;
    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;
  }
  else {
    mm = m12;
    V3 = mm/(2.*m3);
  }

  gdouble ch = MIN(c, 1. - c);
  if (ch < V1)
    alpha = pow (pr*ch, 1./3.);
  else if (ch < V2)
    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;
  else if (ch < V3) {
    gdouble p = 2.*m1*m2;
    gdouble q = 3.*m1*m2*(m12 - 2.*m3*ch)/2.;
    gdouble p12 = sqrt (p);
    gdouble teta = acos(q/(p*p12))/3.;
    gdouble cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;
  }
  else if (m12 < m3)
    alpha = m3*ch + mm/2.;
  else {
    gdouble p = m1*(m2 + m3) + m2*m3 - 1./4.;
    gdouble q = 3.*m1*m2*m3*(1./2. - ch)/2.;
    gdouble p12 = sqrt(p);
    gdouble teta = acos(q/(p*p12))/3.;
    gdouble cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }
  if (c > 1./2.) alpha = 1. - alpha;

  if (m->x < 0.)
    alpha += m->x;
  if (m->y < 0.)
    alpha += m->y;
  if (m->z < 0.)
    alpha += m->z;

  return alpha;
}

/**
 * gfs_plane_center:
 * @m: normal to the plane.
 * @alpha: plane constant.
 * @a: volume of cell fraction.
 * @p: a #FttVector.
 *
 * Fills @p with the position of the center of mass of the fraction of
 * a cubic cell lying under the plane (@m,@alpha).
 */
void gfs_plane_center (const FttVector * m, gdouble alpha, gdouble a, FttVector * p)
{
  FttVector n;
  gdouble b, amax;

  g_return_if_fail (m != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (a >= 0. && a <= 1.);

  if (fabs (m->x) < EPS) {
    FttVector q;
    n.x = m->y;
    n.y = m->z;
    gfs_line_center (&n, alpha, a, &q);
    p->x = 0.5;
    p->y = q.x;
    p->z = q.y;
    return;
  }
  if (fabs (m->y) < EPS) {
    FttVector q;
    n.x = m->z;
    n.y = m->x;
    gfs_line_center (&n, alpha, a, &q);
    p->x = q.y;
    p->y = 0.5;
    p->z = q.x;
    return;
  }
  if (fabs (m->z) < EPS) {
    gfs_line_center (m, alpha, a, p);
    p->z = 0.5;
    return;
  }

  n = *m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }  
  if (n.z < 0.) {
    alpha -= n.z;
    n.z = - n.z;
  }  

  if (alpha <= 0. || a == 0.) {
    p->x = p->y = p->z = 0.;
    return;
  }

  if (alpha >= n.x + n.y + n.z || a == 1.) {
    p->x = p->y = p->z = 0.5;
    return;
  }

  amax = n.x + n.y + n.z;
  p->x = p->y = p->z = alpha*alpha*alpha*alpha;

  b = alpha - n.x;
  if (b > 0.) {
    p->x -= b*b*b*(3.*n.x + alpha);
    p->y -= b*b*b*b;
    p->z -= b*b*b*b;
  }
  b = alpha - n.y;
  if (b > 0.) {
    p->y -= b*b*b*(3.*n.y + alpha);
    p->x -= b*b*b*b;
    p->z -= b*b*b*b;
  }
  b = alpha - n.z;
  if (b > 0.) {
    p->z -= b*b*b*(3.*n.z + alpha);
    p->x -= b*b*b*b;
    p->y -= b*b*b*b;
  }

  amax = alpha - amax;
  b = amax + n.x;
  if (b > 0.) {
    p->y += b*b*b*(3.*n.y + alpha - n.z);
    p->z += b*b*b*(3.*n.z + alpha - n.y);
    p->x += b*b*b*b;
  }
  b = amax + n.y;
  if (b > 0.) {
    p->x += b*b*b*(3.*n.x + alpha - n.z);
    p->z += b*b*b*(3.*n.z + alpha - n.x);
    p->y += b*b*b*b;
  }
  b = amax + n.z;
  if (b > 0.) {
    p->x += b*b*b*(3.*n.x + alpha - n.y);
    p->y += b*b*b*(3.*n.y + alpha - n.x);
    p->z += b*b*b*b;
  }

  b = 24.*n.x*n.y*n.z*a;
  p->x /= b*n.x; p->y /= b*n.y; p->z /= b*n.z;

  if (m->x < 0.) p->x = 1. - p->x;
  if (m->y < 0.) p->y = 1. - p->y;
  if (m->z < 0.) p->z = 1. - p->z;
}

/**
 * gfs_plane_area_center:
 * @m: normal to the plane.
 * @alpha: plane constant.
 * @p: a #FttVector.
 *
 * Fills @p with the position of the center of area of the fraction of
 * a cubic cell lying under the plane (@m,@alpha).
 *
 * Returns: the area of the facet.
 */
gdouble gfs_plane_area_center (const FttVector * m, gdouble alpha, FttVector * p)
{
  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (p != NULL, 0.);

  if (fabs (m->x) < EPS) {
    FttVector n, q;
    n.x = m->y;
    n.y = m->z;
    gdouble area = gfs_line_area_center (&n, alpha, &q);
    p->x = 0.5;
    p->y = q.x;
    p->z = q.y;
    return area;
  }
  if (fabs (m->y) < EPS) {
    FttVector n, q;
    n.x = m->z;
    n.y = m->x;
    gdouble area = gfs_line_area_center (&n, alpha, &q);
    p->x = q.y;
    p->y = 0.5;
    p->z = q.x;
    return area;
  }
  if (fabs (m->z) < EPS) {
    gdouble area = gfs_line_area_center (m, alpha, p);
    p->z = 0.5;
    return area;
  }

  FttVector n = *m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }  
  if (n.z < 0.) {
    alpha -= n.z;
    n.z = - n.z;
  }

  gdouble amax = n.x + n.y + n.z;
  if (alpha <= 0. || alpha >= amax) {
    p->x = p->y = p->z = 0.;
    return 0.;
  }

  gdouble area = alpha*alpha;
  p->x = p->y = p->z = area*alpha;

  gdouble b = alpha - n.x;
  if (b > 0.) {
    area -= b*b;
    p->x -= b*b*(2.*n.x + alpha);
    p->y -= b*b*b;
    p->z -= b*b*b;
  }
  b = alpha - n.y;
  if (b > 0.) {
    area -= b*b;
    p->y -= b*b*(2.*n.y + alpha);
    p->x -= b*b*b;
    p->z -= b*b*b;
  }
  b = alpha - n.z;
  if (b > 0.) {
    area -= b*b;
    p->z -= b*b*(2.*n.z + alpha);
    p->x -= b*b*b;
    p->y -= b*b*b;
  }

  amax = alpha - amax;
  b = amax + n.x;
  if (b > 0.) {
    area += b*b;
    p->y += b*b*(2.*n.y + alpha - n.z);
    p->z += b*b*(2.*n.z + alpha - n.y);
    p->x += b*b*b;
  }
  b = amax + n.y;
  if (b > 0.) {
    area += b*b;
    p->x += b*b*(2.*n.x + alpha - n.z);
    p->z += b*b*(2.*n.z + alpha - n.x);
    p->y += b*b*b;
  }
  b = amax + n.z;
  if (b > 0.) {
    area += b*b;
    p->x += b*b*(2.*n.x + alpha - n.y);
    p->y += b*b*(2.*n.y + alpha - n.x);
    p->z += b*b*b;
  }

  area *= 3.;
  p->x /= area*n.x;
  p->y /= area*n.y;
  p->z /= area*n.z;

  THRESHOLD (p->x);
  THRESHOLD (p->y);
  THRESHOLD (p->z);

  if (m->x < 0.) p->x = 1. - p->x;
  if (m->y < 0.) p->y = 1. - p->y;
  if (m->z < 0.) p->z = 1. - p->z;

  return area*sqrt (1./(n.x*n.x*n.y*n.y) + 1./(n.x*n.x*n.z*n.z) + 1./(n.z*n.z*n.y*n.y))/6.;
}
#endif /* 3D */

/**
 * gfs_youngs_gradient:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 * @g: a #FttVector.
 *
 * Fills @g with the Youngs-averaged gradients of @v 
 * normalised by the size of @cell.
 */
void gfs_youngs_gradient (FttCell * cell, GfsVariable * v, FttVector * g)
{
  static FttDirection d[(FTT_DIMENSION - 1)*4][FTT_DIMENSION] = {
#if FTT_2D
    {FTT_RIGHT, FTT_TOP}, {FTT_LEFT, FTT_TOP}, {FTT_LEFT, FTT_BOTTOM}, {FTT_RIGHT, FTT_BOTTOM}
#else  /* 3D */
    {FTT_RIGHT, FTT_TOP, FTT_FRONT}, {FTT_LEFT, FTT_TOP, FTT_FRONT}, 
    {FTT_LEFT, FTT_BOTTOM, FTT_FRONT}, {FTT_RIGHT, FTT_BOTTOM, FTT_FRONT},
    {FTT_RIGHT, FTT_TOP, FTT_BACK}, {FTT_LEFT, FTT_TOP, FTT_BACK}, 
    {FTT_LEFT, FTT_BOTTOM, FTT_BACK}, {FTT_RIGHT, FTT_BOTTOM, FTT_BACK},
#endif /* 3D */
  };
  gdouble u[(FTT_DIMENSION - 1)*4];
  guint i;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (v != NULL);
  g_return_if_fail (g != NULL);

  for (i = 0; i < (FTT_DIMENSION - 1)*4; i++)
    u[i] = gfs_cell_corner_value (cell, d[i], v, -1);

#if FTT_2D
  g->x = (u[0] + u[3] - u[1] - u[2])/2.;
  g->y = (u[0] + u[1] - u[2] - u[3])/2.;
#else  /* 3D */
  g->x = (u[0] + u[3] + u[4] + u[7] - u[1] - u[2] - u[5] - u[6])/4.;
  g->y = (u[0] + u[1] + u[4] + u[5] - u[2] - u[3] - u[6] - u[7])/4.;
  g->z = (u[0] + u[1] + u[2] + u[3] - u[4] - u[5] - u[6] - u[7])/4.;
#endif /* 3D */
}

/**
 * gfs_vof_plane_interpolate:
 * @cell: a #FttCell containing location @p.
 * @p: the center of the virtual cell.
 * @level: the level of the virtual cell.
 * @t: a #GfsVariableTracerVOF.
 * @m: a #FttVector.
 *
 * Computes the equation @m.x = alpha of the volume fraction plane of
 * a virtual cell at @level centered on @p.
 *
 * Returns: alpha for the virtual cell.
 */
gdouble gfs_vof_plane_interpolate (FttCell * cell,
				   FttVector * p,
				   guint level,
				   GfsVariableTracerVOF * t,
				   FttVector * m)
{
  guint l = ftt_cell_level (cell);

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (l <= level, 0.);
  g_return_val_if_fail (t != NULL, 0.);
  g_return_val_if_fail (m != NULL, 0.);

  GfsVariable * v = GFS_VARIABLE (t);
  gdouble f = GFS_VALUE (cell, v);
  g_return_val_if_fail (!GFS_IS_FULL (f), 0.);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&m->x)[c] = GFS_VALUE (cell, t->m[c]);

  gdouble alpha = GFS_VALUE (cell, t->alpha);
  if (l < level) {
    gdouble h = ftt_level_size (level);
    gdouble H = ftt_cell_size (cell);
    FttVector q;
    
    ftt_cell_pos (cell, &q);
    alpha *= H;
    for (c = 0; c < FTT_DIMENSION; c++)
      alpha -= (&m->x)[c]*((&p->x)[c] - h/2. - (&q.x)[c] + H/2);
    alpha /= h;
  }
  return alpha;
}

/**
 * gfs_vof_interpolate:
 * @cell: a #FttCell containing location @p.
 * @p: the center of the virtual cell.
 * @level: the level of the virtual cell.
 * @t: a #GfsVariableTracerVOF.
 *
 * Computes the volume fraction of a virtual cell at @level centered
 * on @p.
 *
 * Returns: the volume fraction of the virtual cell.
 */
gdouble gfs_vof_interpolate (FttCell * cell,
			     FttVector * p,
			     guint level,
			     GfsVariableTracerVOF * t)
{
  guint l = ftt_cell_level (cell);

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (l <= level, 0.);
  g_return_val_if_fail (t != NULL, 0.);

  GfsVariable * v = GFS_VARIABLE (t);
  gdouble f = GFS_VALUE (cell, v);
  if (l == level || GFS_IS_FULL (f))
    return f;
  else {
    FttVector m;
    gdouble alpha = gfs_vof_plane_interpolate (cell, p, level, t, &m);
    return gfs_plane_volume (&m, alpha);
  }
}

/**
 * Volume-Of-Fluid advection.
 * \beginobject{GfsVariableTracerVOF}
 */

#if FTT_2D
# define F(x,y,z) f[x][y]
#else
# define F(x,y,z) f[x][y][z]
#endif

static void stencil (FttCell * cell, GfsVariable * v, gdouble F(3,3,3))
{
  gdouble h = ftt_cell_size (cell);
  guint level = ftt_cell_level (cell);
  FttVector p;
  gint x, y, z = 0;
  
  F(1,1,1) = GFS_VALUE (cell, v);
  ftt_cell_pos (cell, &p);
#if !FTT_2D
  for (z = -1; z <= 1; z++)
#endif
    for (x = -1; x <= 1; x++)
      for (y = -1; y <= 1; y++)
	if (x != 0 || y != 0 || z != 0) {
	  FttVector o;
	  o.x = p.x + h*x; o.y = p.y + h*y; o.z = p.z + h*z;
	  FttCell * neighbor = gfs_domain_boundary_locate (v->domain, o, level, NULL);
	  if (neighbor)
	    F(x + 1, y + 1, z + 1) =
	      gfs_vof_interpolate (neighbor, &o, level, GFS_VARIABLE_TRACER_VOF (v));
	  else
	    F(x + 1, y + 1, z + 1) = -1.;
	}
  /* boundary conditions (symmetry) */
#if FTT_2D
  for (x = 0; x <= 2; x++) {
    if (f[x][0] < 0.) f[x][0] = f[x][1];
    if (f[x][2] < 0.) f[x][2] = f[x][1];
  }
  for (y = 0; y <= 2; y++) {
    if (f[0][y] < 0.) f[0][y] = f[1][y];
    if (f[2][y] < 0.) f[2][y] = f[1][y];
  }
#else /* 3D */
  for (x = 0; x <= 2; x++)
    for (y = 0; y <= 2; y++) {
      if (f[x][y][0] < 0.) f[x][y][0] = f[x][y][1];
      if (f[x][y][2] < 0.) f[x][y][2] = f[x][y][1];
    }
  for (x = 0; x <= 2; x++)
    for (z = 0; z <= 2; z++) {
      if (f[x][0][z] < 0.) f[x][0][z] = f[x][1][z];
      if (f[x][2][z] < 0.) f[x][2][z] = f[x][1][z];
    }
  for (z = 0; z <= 2; z++)
    for (y = 0; y <= 2; y++) {
      if (f[0][y][z] < 0.) f[0][y][z] = f[1][y][z];
      if (f[2][y][z] < 0.) f[2][y][z] = f[1][y][z];
    }
#endif /* 3D */
}

static void youngs_normal (FttCell * cell, GfsVariable * v, FttVector * n)
{
  gdouble F(3,3,3);

  stencil (cell, v, f);
#if FTT_2D
  n->x = (f[0][2] + 2.*f[0][1] + f[0][0] - 2.*f[2][1] - f[2][2] - f[2][0])/8.;
  n->y = (f[2][0] + 2.*f[1][0] + f[0][0] - 2.*f[1][2] - f[2][2] - f[0][2])/8.;
  n->z = 0.;
#else  /* 3D */
  gdouble mm1 = f[0][0][0] + f[0][0][2] + f[0][2][0] + f[0][2][2] +
    2.*(f[0][0][1] + f[0][2][1] + f[0][1][0] + f[0][1][2]) + 
    4.*f[0][1][1];
  gdouble mm2 = f[2][0][0] + f[2][0][2] + f[2][2][0] + f[2][2][2] + 
    2.*(f[2][0][1] + f[2][2][1] + f[2][1][0] + f[2][1][2]) + 
    4.*f[2][1][1];
  n->x = (mm1 - mm2)/32.;
    
  mm1 = f[0][0][0] + f[0][0][2] + f[2][0][0] + f[2][0][2] + 
    2.*(f[0][0][1] + f[2][0][1] + f[1][0][0] + f[1][0][2]) + 
    4.*f[1][0][1];
  mm2 = f[0][2][0] + f[0][2][2] + f[2][2][0] + f[2][2][2] + 
    2.*(f[0][2][1] + f[2][2][1] + f[1][2][0] + f[1][2][2]) + 
    4.*f[1][2][1];
  n->y = (mm1 - mm2)/32.;
                  
  mm1 = f[0][0][0] + f[0][2][0] + f[2][0][0] + f[2][2][0] +
    2.*(f[0][1][0] + f[2][1][0] + f[1][0][0] + f[1][2][0]) + 
    4.*f[1][1][0];
  mm2 = f[0][0][2] + f[0][2][2] + f[2][0][2] + f[2][2][2] + 
    2.*(f[0][1][2] + f[2][1][2] + f[1][0][2] + f[1][2][2]) + 
    4.*f[1][1][2];
  n->z = (mm1 - mm2)/32.;
#endif /* 3D */
}

#if FTT_2D
# include "myc2d.h"
#else
# include "myc.h"
#endif

static void myc_normal (FttCell * cell, GfsVariable * v, FttVector * n)
{
  gdouble F(3,3,3);  

  stencil (cell, v, f);
  mycs (f, &n->x);
#if FTT_2D
  n->z = 0.;
#endif
}

static void vof_plane (FttCell * cell, GfsVariable * v)
{
  if (FTT_CELL_IS_LEAF (cell)) {
    GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
    gdouble f = GFS_VALUE (cell, v);
    FttComponent c;

    THRESHOLD (f);
    if (GFS_IS_FULL (f)) {
      for (c = 1; c < FTT_DIMENSION; c++)
	GFS_VALUE (cell, t->m[c]) = 0.;
      GFS_VALUE (cell, t->m[0]) = 1.;
      GFS_VALUE (cell, t->alpha) = f;
    }
    else {
      FttVector m;
      gdouble n = 0.;

      myc_normal (cell, v, &m);
      for (c = 0; c < FTT_DIMENSION; c++)
	n += fabs ((&m.x)[c]);
      if (n > 0.)
	for (c = 0; c < FTT_DIMENSION; c++)
	  (&m.x)[c] /= n;
      else /* fixme: this is a small fragment */
	m.x = 1.;
      for (c = 0; c < FTT_DIMENSION; c++)
	GFS_VALUE (cell, t->m[c]) = (&m.x)[c];
      GFS_VALUE (cell, t->alpha) = gfs_plane_alpha (&m, f);
    }
  }
}

static void no_coarse_fine (FttCell * cell,  GfsVariable * v) {}

static void allocate_normal_alpha (GfsVariableTracerVOF * t)
{
  GfsVariable * v = GFS_VARIABLE (t);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++) {
    static gchar index[][2] = {"x", "y", "z"};
    gchar * name = g_strdup_printf ("%s_%s", v->name, index[c]);
    gchar * description = 
      g_strdup_printf ("%s-component of the normal to the interface defined by %s",
		       index[c], v->name);
    t->m[c] = gfs_domain_get_or_add_variable (v->domain, name, description);
    t->m[c]->fine_coarse = t->m[c]->coarse_fine = no_coarse_fine;
    g_free (name);
    g_free (description);
  }
  gchar * name = g_strdup_printf ("%s_alpha", v->name);
  gchar * description = 
    g_strdup_printf ("\"alpha\" for the interface defined by %s", v->name);
  t->alpha = gfs_domain_get_or_add_variable (v->domain, name, description);
  t->alpha->fine_coarse = t->alpha->coarse_fine = no_coarse_fine;
  g_free (name);
  g_free (description);
}

static void variable_tracer_vof_update (GfsVariable * v, GfsDomain * domain)
{
  GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
  gfs_domain_cell_traverse (domain,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) v->fine_coarse, v);
  gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, v);
  GSList * j = t->concentrations->items;
  while (j) {
    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, j->data);
    j = j->next;
  }
    
  /* update normals and alpha */
  guint l, depth = gfs_domain_depth (domain);
  FttComponent c;
  for (l = 0; l <= depth; l++) {
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, l,
			      (FttCellTraverseFunc) vof_plane, v);
    for (c = 0; c < FTT_DIMENSION; c++)
      gfs_domain_bc (domain, FTT_TRAVERSE_LEVEL, l, t->m[c]);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEVEL, l, t->alpha);
  }
}

static gboolean variable_tracer_vof_event (GfsEvent * event, 
					   GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_tracer_vof_class ())->parent_class)->event)
      (event, sim)) {
    (* GFS_VARIABLE_TRACER_VOF_CLASS (GTS_OBJECT (event)->klass)->update) (GFS_VARIABLE (event),
									   GFS_DOMAIN (sim));
    return TRUE;
  }
  return FALSE;
}

static void variable_tracer_vof_destroy (GtsObject * o)
{
  GfsVariableTracerVOF * v = GFS_VARIABLE_TRACER_VOF (o);

  if (v->alpha) {
    FttComponent c;
    for (c = 0; c < FTT_DIMENSION; c++)
      gts_object_destroy (GTS_OBJECT (v->m[c]));
    gts_object_destroy (GTS_OBJECT (v->alpha));
  }
  gts_object_destroy (GTS_OBJECT (v->concentrations));

  (* GTS_OBJECT_CLASS (gfs_variable_tracer_vof_class ())->parent_class->destroy) (o);
}

static void variable_tracer_vof_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_tracer_vof_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (GFS_VARIABLE_TRACER (*o)->advection.cfl > 0.5) {
    gts_file_error (fp, "cfl `%g' is out of range `]0,0.5]'", 
		    GFS_VARIABLE_TRACER (*o)->advection.cfl);
    return;
  }  

  allocate_normal_alpha (GFS_VARIABLE_TRACER_VOF (*o));
}

static void variable_tracer_vof_class_init (GtsObjectClass * klass)
{
  GFS_VARIABLE_TRACER_VOF_CLASS (klass)->update = variable_tracer_vof_update;
  GFS_EVENT_CLASS (klass)->event = variable_tracer_vof_event;
  klass->destroy = variable_tracer_vof_destroy;
  klass->read = variable_tracer_vof_read;
}

static void vof_coarse_fine (FttCell * parent, GfsVariable * v)
{
  GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
  gdouble f = GFS_VALUE (parent, v);
  FttCellChildren child;
  guint i;
  
  ftt_cell_children (parent, &child);
  if (GFS_IS_FULL (f)) {
    for (i = 0; i < FTT_CELLS; i++) 
      if (child.c[i]) {
	FttComponent c;
	GFS_VALUE (child.c[i], v) = f;
	for (c = 1; c < FTT_DIMENSION; c++)
	  GFS_VALUE (child.c[i], t->m[c]) = 0.;
	GFS_VALUE (child.c[i], t->m[0]) = 1.;
	GFS_VALUE (child.c[i], t->alpha) = f;
      }
    
    /* second-order interpolation for concentrations in full cells */
    GSList * j = t->concentrations->items;
    while (j) {
      gfs_cell_coarse_fine (parent, j->data);
      j = j->next;
    }
  }
  else {
    /* interfacial cell */
    gdouble alpha = GFS_VALUE (parent, t->alpha);
    FttVector m;
    
    for (i = 0; i < FTT_DIMENSION; i++)
      (&m.x)[i] = GFS_VALUE (parent, t->m[i]);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i]) {
	/* fixme: this is not mass conserving */
	gdouble alpha1 = alpha;
	FttComponent c;
	FttVector p;
	ftt_cell_relative_pos (child.c[i], &p);
	for (c = 0; c < FTT_DIMENSION; c++) {
	  alpha1 -= (&m.x)[c]*(0.25 + (&p.x)[c]);
	  GFS_VALUE (child.c[i], t->m[c]) = (&m.x)[c];
	}
	GFS_VALUE (child.c[i], v) = gfs_plane_volume (&m, 2.*alpha1);
	GFS_VALUE (child.c[i], t->alpha) = 2.*alpha1;
      }

    /* firs-order interpolation for concentrations in interfacial cells */
    GSList * j = t->concentrations->items;
    while (j) {
      GfsVariable * v = j->data;
      for (i = 0; i < FTT_CELLS; i++)
	if (child.c[i])
	  GFS_VALUE (child.c[i], v) = GFS_VALUE (parent, v);
      j = j->next;
    }
  }
}

static void vof_fine_coarse (FttCell * parent, GfsVariable * v)
{
  gdouble val = 0., sa = 0.;
  guint i;
  FttCellChildren child;

  ftt_cell_children (parent, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i] && GFS_HAS_DATA (child.c[i], v)) {
      val += GFS_VALUE (child.c[i], v);
      sa += 1.;
    }
  if (sa > 0.)
    GFS_VALUE (parent, v) = val/sa;
  else
    GFS_VALUE (parent, v) = GFS_NODATA;

  GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
  gdouble f = GFS_VALUE (parent, v);
  FttComponent c;

  if (GFS_IS_FULL (f)) {
    for (c = 1; c < FTT_DIMENSION; c++)
      GFS_VALUE (parent, t->m[c]) = 0.;
    GFS_VALUE (parent, t->m[0]) = 1.;
    GFS_VALUE (parent, t->alpha) = f;
  }
  else {
    FttVector m = {0., 0., 0.};
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i]) {
	gdouble f = GFS_VALUE (child.c[i], v);
	gdouble a = f*(1. - f);

	for (c = 0; c < FTT_DIMENSION; c++)
	  (&m.x)[c] += a*GFS_VALUE (child.c[i], t->m[c]);
      }
    
    gdouble n = 0.;
    for (c = 0; c < FTT_DIMENSION; c++)
      n += fabs ((&m.x)[c]);
    if (n > 0.)
      for (c = 0; c < FTT_DIMENSION; c++)
	(&m.x)[c] /= n;
    else /* fixme: this is a small fragment */
      m.x = 1.;
    for (c = 0; c < FTT_DIMENSION; c++)
      GFS_VALUE (parent, t->m[c]) = (&m.x)[c];
    GFS_VALUE (parent, t->alpha) = gfs_plane_alpha (&m, f);
  }

  GSList * j = t->concentrations->items;
  while (j) {
    GfsVariable * p = j->data;
    gdouble sp = 0.;
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i] && GFS_HAS_DATA (child.c[i], v))
	sp += GFS_VALUE (child.c[i], v)*GFS_VALUE (child.c[i], p);
    if (val > 0.)
      GFS_VALUE (parent, p) = sp/val;
    else
      GFS_VALUE (parent, p) = 0.;
    j = j->next;
  }
}

static void variable_tracer_vof_init (GfsVariable * v)
{
  GFS_EVENT (v)->start = -1;
  GFS_EVENT (v)->istep = G_MAXINT/2;
  v->coarse_fine = vof_coarse_fine;
  v->fine_coarse = vof_fine_coarse;
  //  v->face_value = gfs_vof_face_value;
  GFS_VARIABLE_TRACER (v)->advection.cfl = 0.5;
  GFS_VARIABLE_TRACER_VOF (v)->concentrations = 
    GTS_SLIST_CONTAINER (gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ())));
}

GfsVariableTracerVOFClass * gfs_variable_tracer_vof_class (void)
{
  static GfsVariableTracerVOFClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsVariableTracerVOF",
      sizeof (GfsVariableTracerVOF),
      sizeof (GfsVariableTracerVOFClass),
      (GtsObjectClassInitFunc) variable_tracer_vof_class_init,
      (GtsObjectInitFunc) variable_tracer_vof_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_tracer_class ()), &info);
  }

  return klass;
}

typedef struct {
  GfsAdvectionParams * par, vpar;
  GfsVariable * u, * du[FTT_DIMENSION - 1], * vof;
  FttComponent c;
  GfsDomain * domain;
  GfsFunction * sink;
  guint depth, too_coarse;
} VofParms;

static gdouble plane_volume_shifted (FttVector m, gdouble alpha, FttVector p[2])
{
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++) {
    alpha -= (&m.x)[c]*(&p[0].x)[c];
    (&m.x)[c] *= (&p[1].x)[c] - (&p[0].x)[c];
  }
  return gfs_plane_volume (&m, alpha);
}

static gdouble fine_fraction (FttCellFace * face, GfsVariable * v, gdouble un, FttVector q[2])
{
  gdouble f = GFS_VALUE (face->cell, v);
  if (f == 0. || f == 1.)
    return f;
  else {
    FttComponent c;
    FttVector m, q1[2];
    gdouble alpha = GFS_VALUE (face->cell, GFS_VARIABLE_TRACER_VOF (v)->alpha);

    for (c = 0; c < FTT_DIMENSION; c++)
      (&m.x)[c] = GFS_VALUE (face->cell, GFS_VARIABLE_TRACER_VOF (v)->m[c]);
    c = face->d/2;
    if (face->d % 2 != 0) {
      (&m.x)[c] = - (&m.x)[c];
      alpha += (&m.x)[c];
    }

    q1[0] = q[0]; q1[1] = q[1];
    (&q1[0].x)[c] = 1. - un;
    return plane_volume_shifted (m, alpha, q1);
  }
}

static gdouble coarse_fraction (FttCellFace * face, GfsVariable * v, gdouble un, FttVector q[2])
{
  gdouble f = GFS_VALUE (face->neighbor, v);
  if (f == 0. || f == 1.)
    return f;
  else {
    FttComponent c;
    FttVector m, o, q1[2];
    gdouble alpha = GFS_VALUE (face->neighbor, GFS_VARIABLE_TRACER_VOF (v)->alpha);
    
    for (c = 0; c < FTT_DIMENSION; c++)
      (&m.x)[c] = GFS_VALUE (face->neighbor, GFS_VARIABLE_TRACER_VOF (v)->m[c]);
    
    if (!FTT_FACE_DIRECT (face)) {
      (&m.x)[face->d/2] = - (&m.x)[face->d/2];
      alpha += (&m.x)[face->d/2];
    }
    
    /* shift interface perpendicularly */
    q1[0] = q[0]; q1[1] = q[1];
    ftt_cell_relative_pos (face->cell, &o);
    for (c = 0; c < FTT_DIMENSION; c++)
      if (c != face->d/2) {
	(&q1[0].x)[c] = (&o.x)[c] + 1./4. + (&q1[0].x)[c]/2.;
	(&q1[1].x)[c] = (&o.x)[c] + 1./4. + (&q1[1].x)[c]/2.;
      }
    (&q1[1].x)[face->d/2] = un;
    return plane_volume_shifted (m, alpha, q1);
  }
}

#define TOO_COARSE(cell) (GFS_VALUE (cell, p->par->fv))

/* Marks coarse cells which should be refined because an interface in
   a neighboring finer cell will be advected into them */
static void face_too_coarse (FttCellFace * face, VofParms * p)
{
  if (ftt_face_type (face) == FTT_FINE_COARSE) {
    gdouble un = GFS_FACE_NORMAL_VELOCITY (face);
    if (!FTT_FACE_DIRECT (face))
      un = - un;
    if (un > 0.) {
      gdouble f = GFS_VALUE (face->neighbor, p->par->v);
      if (GFS_IS_FULL (f)) {
	FttVector q[2] = {{0., 0., 0.},{1., 1., 1.}};
	if (fine_fraction (face, p->vof, un*p->par->dt/ftt_cell_size (face->cell), q) != f) {
	  p->too_coarse++;
	  TOO_COARSE (face->neighbor) = TRUE;
	}
      }
    }
  }
}

static void vof_cell_fine_init (FttCell * parent, VofParms * p)
{
  gfs_cell_fine_init (parent, p->domain);

  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    FttDirection od = FTT_OPPOSITE_DIRECTION (d);
    FttCellChildren dchild;
    guint i, n = ftt_cell_children_direction (parent, d, &dchild);
    for (i = 0; i < n; i++)
      if (dchild.c[i]) {
	FttCell * neighbor = ftt_cell_neighbor (dchild.c[i], d);
	if (neighbor)
	  GFS_STATE (dchild.c[i])->f[d].un = GFS_STATE (neighbor)->f[od].un;
      }
  }

  FttCellChildren child;
  gdouble div[FTT_CELLS], P[FTT_CELLS];
  guint n;
  ftt_cell_children (parent, &child);
  for (n = 0; n < FTT_CELLS; n++) {
    div[n] = 0.;
    if (child.c[n]) {
      GFS_VALUE (child.c[n], p->vpar.v) = GFS_VALUE (parent, p->vpar.v);
      FttComponent c;
      for (c = 0; c < FTT_DIMENSION; c++)
	div[n] += GFS_STATE (child.c[n])->f[2*c].un - GFS_STATE (child.c[n])->f[2*c + 1].un;
    }
  }

#if FTT_2D
  P[0] = 0.;
  P[1] = (3.*div[1] + div[2])/4. + div[3]/2.;
  P[2] = (div[1] + 3.*div[2])/4. + div[3]/2.;
  P[3] = (div[1] + div[2])/2. + div[3];
  if (child.c[0]) {
    GFS_STATE (child.c[0])->f[0].un = P[1] - P[0];
    GFS_STATE (child.c[0])->f[3].un = P[0] - P[2];
  }
  if (child.c[1]) {
    GFS_STATE (child.c[1])->f[1].un = P[1] - P[0];
    GFS_STATE (child.c[1])->f[3].un = P[1] - P[3];
  }
  if (child.c[2]) {
    GFS_STATE (child.c[2])->f[0].un = P[3] - P[2];
    GFS_STATE (child.c[2])->f[2].un = P[0] - P[2];
  }
  if (child.c[3]) {
    GFS_STATE (child.c[3])->f[1].un = P[3] - P[2];
    GFS_STATE (child.c[3])->f[2].un = P[1] - P[3];
  }
#else /* 3D */
  static gdouble m[7][7] = {{7./12.,5./24.,3./8.,5./24.,3./8.,1./4.,1./3.},
			    {5./24.,7./12.,3./8.,5./24.,1./4.,3./8.,1./3.}, 
			    {3./8.,3./8.,3./4.,1./4.,3./8.,3./8.,1./2.}, 
			    {5./24.,5./24.,1./4.,7./12.,3./8.,3./8.,1./3.}, 
			    {3./8.,1./4.,3./8.,3./8.,3./4.,3./8.,1./2.}, 
			    {1./4.,3./8.,3./8.,3./8.,3./8.,3./4.,1./2.}, 
			    {1./3.,1./3.,1./2.,1./3.,1./2.,1./2.,5./6.}};
  P[0] = 0.;
  guint i, j;
  for (i = 0; i < 7; i++) {
    P[i + 1] = 0.;
    for (j = 0; j < 7; j++)
      P[i + 1] += m[i][j]*div[j + 1];
  }
  if (child.c[0]) {
    GFS_STATE (child.c[0])->f[0].un = P[1] - P[0];
    GFS_STATE (child.c[0])->f[3].un = P[0] - P[2];
    GFS_STATE (child.c[0])->f[5].un = P[0] - P[4];
  }
  if (child.c[1]) {
    GFS_STATE (child.c[1])->f[1].un = P[1] - P[0];
    GFS_STATE (child.c[1])->f[3].un = P[1] - P[3];
    GFS_STATE (child.c[1])->f[5].un = P[1] - P[5];
  }
  if (child.c[2]) {
    GFS_STATE (child.c[2])->f[0].un = P[3] - P[2];
    GFS_STATE (child.c[2])->f[2].un = P[0] - P[2];
    GFS_STATE (child.c[2])->f[5].un = P[2] - P[6];
  }
  if (child.c[3]) {
    GFS_STATE (child.c[3])->f[1].un = P[3] - P[2];
    GFS_STATE (child.c[3])->f[2].un = P[1] - P[3];
    GFS_STATE (child.c[3])->f[5].un = P[3] - P[7];
  }
  if (child.c[4]) {
    GFS_STATE (child.c[4])->f[0].un = P[5] - P[4];
    GFS_STATE (child.c[4])->f[3].un = P[4] - P[6];
    GFS_STATE (child.c[4])->f[4].un = P[0] - P[4];
  }
  if (child.c[5]) {
    GFS_STATE (child.c[5])->f[1].un = P[5] - P[4];
    GFS_STATE (child.c[5])->f[3].un = P[5] - P[7];
    GFS_STATE (child.c[5])->f[4].un = P[1] - P[5];
  }
  if (child.c[6]) {
    GFS_STATE (child.c[6])->f[0].un = P[7] - P[6];
    GFS_STATE (child.c[6])->f[2].un = P[4] - P[6];
    GFS_STATE (child.c[6])->f[4].un = P[2] - P[6];
  }
  if (child.c[7]) {
    GFS_STATE (child.c[7])->f[1].un = P[7] - P[6];
    GFS_STATE (child.c[7])->f[2].un = P[5] - P[7];
    GFS_STATE (child.c[7])->f[4].un = P[3] - P[7];
  }
#endif /* 3D */
}

/* Same as vof_cell_fine_init() but initialisation of MAC velocities
   is done through a GfsVariableStreamFunction */
static void vof_cell_fine_init_with_streamfunction (FttCell * parent, VofParms * p)
{
  gfs_cell_fine_init (parent, p->domain);

  FttCellChildren child;
  guint n;
  ftt_cell_children (parent, &child);
  for (n = 0; n < FTT_CELLS; n++) {
    g_assert (child.c[n]);
    GFS_VALUE (child.c[n], p->vpar.v) = GFS_VALUE (parent, p->vpar.v);
  }
}

static void refine_too_coarse (FttCell * cell, VofParms * p)
{
  if (TOO_COARSE (cell)) {
    guint level = ftt_cell_level (cell);

    TOO_COARSE (cell) = FALSE;
    ftt_cell_refine_corners (cell, (FttCellInitFunc) p->domain->cell_init, p);
    ftt_cell_refine_single (cell, (FttCellInitFunc) p->domain->cell_init, p);
    if (level + 1 > p->depth)
      p->depth = level + 1;
  }
}

/* refine cells which would lead to a loss of resolution at the interface */
static void fix_too_coarse (GfsDomain * domain, VofParms * p)
{
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, p->par->fv);
  p->depth = 0;
  p->domain = domain;
  p->too_coarse = 0;
  gfs_domain_face_traverse (domain, p->c,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) face_too_coarse, p);
  gboolean streamfunction = FALSE;
#if FTT_2D
  GSList * i = domain->variables;
  while (i && !streamfunction) {
    streamfunction = (GFS_IS_VARIABLE_STREAM_FUNCTION (i->data) != NULL);
    i = i->next;
  }
#endif /* 2D */ 
  domain->cell_init = (FttCellInitFunc) (streamfunction ? 
					 vof_cell_fine_init_with_streamfunction : 
					 vof_cell_fine_init);
  domain->cell_init_data = p;
  if (p->too_coarse > 0)
    gfs_domain_cell_traverse (domain,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) refine_too_coarse, p);
  gfs_all_reduce (domain, p->too_coarse, MPI_UNSIGNED, MPI_SUM);
  if (p->too_coarse > 0)
    gfs_domain_reshape (domain, p->depth);
  domain->cell_init = (FttCellInitFunc) gfs_cell_fine_init;
  domain->cell_init_data = domain;
}

static void concentration_face_values (FttCell * cell, VofParms * p)
{
  GfsStateVector * s = GFS_STATE (cell);
  gdouble size = ftt_cell_size (cell);
  if (p->domain->scale_metric)
      size *= (* p->domain->scale_metric) (p->domain, cell, p->c);
  gdouble unorm = p->par->dt*(s->f[2*p->c].un + s->f[2*p->c + 1].un)/(2.*size);
  if (p->sink)
    unorm += p->par->dt*gfs_function_value (p->sink, cell)/size;
  gdouble g = (* p->par->gradient) (cell, p->c, p->par->v->i);
  gdouble v = GFS_VALUE (cell, p->par->v);
  s->f[2*p->c].v     = v + MIN ((  1. - unorm)/2.,  0.5)*g;
  s->f[2*p->c + 1].v = v + MAX ((- 1. - unorm)/2., -0.5)*g;
  GFS_VALUE (cell, p->par->fv) = 0.;
}

static void vof_flux (FttCellFace * face, VofParms * p)
{
  gdouble size = ftt_cell_size (face->cell);
  gdouble un = GFS_FACE_NORMAL_VELOCITY (face)*p->par->dt/size, dun[FTT_DIMENSION - 1];
  if (p->sink)
    un += gfs_function_face_value (p->sink, face)*p->par->dt/size;
  FttComponent c;

  int n; /* loop over n "horizontal bands" */
  if (GFS_IS_FULL (GFS_VALUE (face->cell, p->vof)) && 
      GFS_IS_FULL (GFS_VALUE (face->neighbor, p->vof))) {
    n = 1; /* non-interfacial cells: use only one band */
    for (c = 0; c < FTT_DIMENSION - 1; c++)
      dun[c] = 0.;
  }
  else {
    for (c = 0; c < FTT_DIMENSION - 1; c++)
      dun[c] = gfs_face_interpolated_value (face, p->du[c]->i)*p->par->dt;
    n = 4; /* interfacial cells: use 4 bands */
  }

  if (!FTT_FACE_DIRECT (face)) {
    un = - un;
    for (c = 0; c < FTT_DIMENSION - 1; c++)
      dun[c] = - dun[c];
  }
  if (fabs (un) > 0.51) {
    FttVector p;
    ftt_face_pos (face, &p);
    g_warning ("CFL (%g) at (%g,%g,%g) is larger than 0.51!", un, p.x, p.y, p.z);
  }

  FttVector q[2] = {{0., 0., 0.},{1., 1., 1.}};
  gdouble flux = 0., unj = un;
  FttComponent ci = FTT_ORTHOGONAL_COMPONENT (p->c);

#if FTT_2D
  gdouble f = gfs_domain_face_fraction (p->vof->domain, face)/n;
#else /* 3D */
  gdouble f = gfs_domain_face_fraction (p->vof->domain, face)/(n*n);
  FttComponent cj = FTT_ORTHOGONAL_COMPONENT (ci);
  int j;
  for (j = 0; j < n; j++) {
    (&q[0].x)[cj] = j/(gdouble) n;        /* front of the band */
    (&q[1].x)[cj] = (j + 1)/(gdouble) n;  /* back of the band */
    /* linear interpolation of the velocity at the center of the band */
    unj = un + (1 - n + 2*j)*dun[1]/(gdouble) (2*n);
#endif /* 3D */

  int i;
  for (i = 0; i < n; i++) {
    (&q[0].x)[ci] = i/(gdouble) n;        /* bottom of the band */
    (&q[1].x)[ci] = (i + 1)/(gdouble) n;  /* top of the band */
    /* linear interpolation of the velocity at the center of the band */
    gdouble uni = unj + (1 - n + 2*i)*dun[0]/(gdouble) (2*n);
    
    switch (ftt_face_type (face)) {
    case FTT_FINE_FINE: {
      if (uni < 0.) {
	FttCell * tmp = face->cell;
	face->cell = face->neighbor;
	face->neighbor = tmp;
	face->d = FTT_OPPOSITE_DIRECTION (face->d);
	uni = - uni;
	unj = - unj;
	un = - un;
	for (c = 0; c < FTT_DIMENSION - 1; c++)
	  dun[c] = - dun[c];
      }
      flux = fine_fraction (face, p->vof, uni, q)*uni*f;
      if (p->par->v == p->vof)
	GFS_VALUE (face->neighbor, p->vpar.fv) += uni*f;
      else
	flux *= GFS_STATE (face->cell)->f[face->d].v;
      GFS_VALUE (face->neighbor, p->par->fv) += flux;	
      break;
    }
    case FTT_FINE_COARSE: {
      flux = uni > 0. ? 
	fine_fraction (face, p->vof, uni, q)*uni*f : 
	coarse_fraction (face, p->vof, -uni/2., q)*uni*f;
      if (p->par->v == p->vof)
	GFS_VALUE (face->neighbor, p->vpar.fv) += uni*f/FTT_CELLS;
      else
	flux *= uni > 0. ?
	  GFS_STATE (face->cell)->f[face->d].v :
	  GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v;
      GFS_VALUE (face->neighbor, p->par->fv) += flux/FTT_CELLS;
      break;
    }
    default:
      g_assert_not_reached ();
    }
    if (p->par->v == p->vof)
      GFS_VALUE (face->cell, p->vpar.fv) -= uni*f;
    GFS_VALUE (face->cell, p->par->fv) -= flux;
  }

#if !FTT_2D
  }
#endif
}

static void initialize_dV (FttCell * cell, GfsVariable * dV)
{
  GFS_VALUE (cell, dV) = 1.;
}

static void reset_fluxes (FttCellFace * face, VofParms * p)
{
  GFS_VALUE (face->cell, p->par->fv) = GFS_VALUE (face->neighbor, p->par->fv) = 0.;
  GFS_VALUE (face->cell, p->vpar.fv) = GFS_VALUE (face->neighbor, p->vpar.fv) = 0.;
}

static void grad_u (FttCell * cell, VofParms * p)
{
  FttComponent c, d = FTT_ORTHOGONAL_COMPONENT (p->c);
  for (c = 0; c < FTT_DIMENSION - 1; c++) {
    GFS_VALUE (cell, p->du[c]) = gfs_center_gradient (cell, d, p->u->i)/ftt_cell_size (cell);
    d = FTT_ORTHOGONAL_COMPONENT (d);
  }
}

static void f_times_dV (FttCell * cell, VofParms * p)
{
  GFS_VALUE (cell, p->par->v) *= GFS_VALUE (cell, p->vpar.v);
}

static void concentration_times_dV (FttCell * cell, VofParms * p)
{
  GFS_VALUE (cell, p->par->v) *= GFS_VALUE (cell, p->vpar.v)*GFS_VALUE (cell, p->vof);
}

static void f_over_dV (FttCell * cell, VofParms * p)
{
  g_assert (GFS_VALUE (cell, p->vpar.v) > 0.);
  gdouble f = GFS_VALUE (cell, p->par->v)/GFS_VALUE (cell, p->vpar.v);
  GFS_VALUE (cell, p->par->v) = f < 1e-10 ? 0. : f > 1. - 1e-10 ? 1. : f;
}

static void concentration_over_dV (FttCell * cell, VofParms * p)
{
  gdouble f = GFS_VALUE (cell, p->vof);
  if (f > 0.)
    GFS_VALUE (cell, p->par->v) = GFS_VALUE (cell, p->par->v)/(GFS_VALUE (cell, p->vpar.v)*f);
  else
    GFS_VALUE (cell, p->par->v) = GFS_NODATA;
}

static void per_vof_volume (FttCell * cell, GfsVariable * v)
{
  GfsVariable * vof = GFS_VARIABLE (GFS_VARIABLE_VOF_CONCENTRATION (v)->vof);
  gdouble f = GFS_VALUE (cell, vof);
  GFS_VALUE (cell, v) = f > 0. ? GFS_VALUE (cell, v)/f : GFS_NODATA;
}

static void per_cell_volume (FttCell * cell, GfsVariable * v)
{
  GfsVariable * vof = GFS_VARIABLE (GFS_VARIABLE_VOF_CONCENTRATION (v)->vof);
  GFS_VALUE (cell, v) *= GFS_VALUE (cell, vof);
}

static void add_sink_velocity (FttCell * cell, VofParms * p)
{
  GFS_VALUE (cell, p->u) += gfs_function_value (p->sink, cell);
}

static void remove_sink_velocity (FttCell * cell, VofParms * p)
{
  GFS_VALUE (cell, p->u) -= gfs_function_value (p->sink, cell);
}

/**
 * gfs_tracer_vof_advection:
 * @domain: a #GfsDomain.
 * @par: the advection parameters.
 *
 * Advects the @v field of @par using the current face-centered (MAC)
 * velocity field.
 */
void gfs_tracer_vof_advection (GfsDomain * domain,
			       GfsAdvectionParams * par)
{
  VofParms p;
  static FttComponent cstart = 0;
  FttComponent c, d;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);
  g_return_if_fail (GFS_IS_VARIABLE_TRACER_VOF (par->v));
  g_return_if_fail (par->cfl <= 0.5);

  gfs_domain_timer_start (domain, "tracer_vof_advection");

  p.par = par;
  p.vof = par->v;
  p.sink = NULL;
  gfs_advection_params_init (&p.vpar);
  for (d = 0; d < FTT_DIMENSION - 1; d++)
    p.du[d] = gfs_temporary_variable (domain);
  p.vpar.v = gfs_temporary_variable (domain);
  p.vpar.fv = gfs_temporary_variable (domain);
  p.vpar.average = par->average;
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) initialize_dV, p.vpar.v);
  par->fv = gfs_temporary_variable (domain);
  GSList * concentrations = GFS_VARIABLE_TRACER_VOF (p.vof)->concentrations->items, * j;
  j = concentrations;
  while (j) {
    GFS_VARIABLE_TRACER (j->data)->advection.fv = gfs_temporary_variable (domain);
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) per_vof_volume, j->data);
    j = j->next;
  }
  for (c = 0; c < FTT_DIMENSION; c++) {
    p.c = (cstart + c) % FTT_DIMENSION;
    fix_too_coarse (domain, &p);
    p.u = gfs_domain_velocity (domain)[p.c];
    gfs_domain_face_traverse (domain, p.c,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) reset_fluxes, &p);
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) grad_u, &p);
    for (d = 0; d < FTT_DIMENSION - 1; d++)
      gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, p.du[d]);
    gfs_domain_face_traverse (domain, p.c,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) vof_flux, &p);
    j = concentrations;
    while (j) {
      GfsAdvectionParams * par = &GFS_VARIABLE_TRACER (j->data)->advection;
      GfsVariable * fv = p.par->fv;
      p.par->v = j->data;
      p.par->fv = par->fv;
      p.par->gradient = par->gradient;
      if (par->sink[0]) {
	p.sink = par->sink[p.c];
	gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) add_sink_velocity, &p);
	gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) grad_u, &p);
	for (d = 0; d < FTT_DIMENSION - 1; d++)
	  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, p.du[d]);
      }
      gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) concentration_face_values, &p);
      gfs_domain_face_bc (domain, p.c, p.par->v);
      gfs_domain_face_traverse (domain, p.c,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttFaceTraverseFunc) vof_flux, &p);
      if (p.sink) {
	gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) remove_sink_velocity, &p);
	p.sink = NULL;
	gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) grad_u, &p);
	for (d = 0; d < FTT_DIMENSION - 1; d++)
	  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, p.du[d]);
      }
      gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) concentration_times_dV, &p);
      gfs_domain_traverse_merged (domain, (GfsMergedTraverseFunc) par->update, par);
      p.par->fv = fv;
      p.par->v = p.vof;
      j = j->next;
    }
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) f_times_dV, &p);
    gfs_domain_traverse_merged (domain, (GfsMergedTraverseFunc) par->update, par);
    gfs_domain_traverse_merged (domain, (GfsMergedTraverseFunc) par->update, &p.vpar);
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) f_over_dV, &p);
    j = concentrations;
    while (j) {
      p.par->v = j->data;
      gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) concentration_over_dV, &p);
      p.par->v = p.vof;
      j = j->next;
    }

    /* update VOF data (normals etc...) */
    (* GFS_VARIABLE_TRACER_VOF_CLASS (GTS_OBJECT (p.par->v)->klass)->update) (p.par->v, domain);
  }
  cstart = (cstart + 1) % FTT_DIMENSION;
  gts_object_destroy (GTS_OBJECT (par->fv));
  par->fv = NULL;
  j = concentrations;
  while (j) {
    gts_object_destroy (GTS_OBJECT (GFS_VARIABLE_TRACER (j->data)->advection.fv));
    GFS_VARIABLE_TRACER (j->data)->advection.fv = NULL;
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) per_cell_volume, j->data);
    j = j->next;
  }
  gts_object_destroy (GTS_OBJECT (p.vpar.v));
  gts_object_destroy (GTS_OBJECT (p.vpar.fv));
  for (d = 0; d < FTT_DIMENSION - 1; d++)
    gts_object_destroy (GTS_OBJECT (p.du[d]));

  gfs_domain_timer_stop (domain, "tracer_vof_advection");
}

static gdouble face_value (FttCell * cell, FttDirection d, GfsVariable * v)
{
  gdouble f = GFS_VALUE (cell, v);

  if (GFS_IS_FULL (f))
    return f;
  else {
    GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
    gdouble alpha = GFS_VALUE (cell, t->alpha);
    FttComponent c;
    FttVector m;
    
    for (c = 0; c < FTT_DIMENSION; c++)
      (&m.x)[c] = GFS_VALUE (cell, t->m[c]);
    (&m.x)[d/2] /= 2.;
    if (d % 2 == 0) {
      (&m.x)[d/2] = -(&m.x)[d/2];
      alpha += (&m.x)[d/2];
    }
    (&m.x)[d/2] /= 2.;
    return gfs_plane_volume (&m, alpha);
  }
}

/**
 * gfs_vof_face_value:
 * @face: a #FttCellFace.
 * @t: a #GfsVariableTracerVOF.
 *
 * Returns: the value of the VOF fraction defined by @t, interpolated
 * on @face.
 */
gdouble gfs_vof_face_value (const FttCellFace * face, GfsVariableTracerVOF * t)
{
  g_return_val_if_fail (face != NULL, 0.);
  g_return_val_if_fail (t != NULL, 0.);

  GfsVariable * v = GFS_VARIABLE (t);
  gdouble vright, vleft = GFS_VALUE (face->cell, v); //face_value (face->cell, face->d, v);
  if (ftt_face_type (face) == FTT_FINE_COARSE) {
    gdouble f = GFS_VALUE (face->neighbor, v);

    if (GFS_IS_FULL (f))
      vright = f;
    else {
      gdouble alpha = GFS_VALUE (face->neighbor, t->alpha);
      FttComponent c;
      FttVector m;

      for (c = 0; c < FTT_DIMENSION; c++)
	(&m.x)[c] = GFS_VALUE (face->neighbor, t->m[c]);

      FttVector p, o;
      ftt_face_pos (face, &p);
      ftt_cell_pos (face->neighbor, &o);
      gdouble h = ftt_cell_size (face->neighbor);

      (&p.x)[face->d/2] += face->d % 2 ? -h/4. : h/4.;
      for (c = 0; c < FTT_DIMENSION; c++)
	alpha -= (&m.x)[c]*(0.25 - ((&p.x)[c] - (&o.x)[c])/h);
      alpha *= 2.;
      //      if (face->d % 2 == 0) {
      //	(&m.x)[face->d/2] = -(&m.x)[face->d/2];
      //	alpha += (&m.x)[face->d/2];
      //      }
      //      (&m.x)[face->d/2] /= 2.;
      vright = gfs_plane_volume (&m, alpha);
#if 0
      if (vright > 0.2 && vright < 0.8) {
	fprintf (stderr, "%d (%g,%g) (%g,%g) %g\n", face->d, p.x, p.y, o.x, o.y, vright);
	g_assert_not_reached ();
      }
#endif
    }
  }
  else
    vright = GFS_VALUE (face->neighbor, v); //face_value (face->neighbor, FTT_OPPOSITE_DIRECTION (face->d), v);
  return (vright + vleft)/2.;
}

/**
 * gfs_vof_facet:
 * @cell: a #FttCell.
 * @t: a #GfsVariableTracerVOF.
 * @p: a #FttVector array (of size 2 in 2D and 6 in 3D)
 * @m: a #FttVector.
 *
 * Fills @p with the coordinates of points defining the
 * VOF-reconstructed interface facet defined by @t.
 *
 * Fills @m with the normal to the interface.
 *
 * Returns: the number of points defining the facet.
 */
guint gfs_vof_facet (FttCell * cell,
		     GfsVariableTracerVOF * t,
		     FttVector * p,
		     FttVector * m)
{
  g_return_val_if_fail (cell != NULL, 0);
  g_return_val_if_fail (t != NULL, 0);
  g_return_val_if_fail (p != NULL, 0);
  g_return_val_if_fail (m != NULL, 0);

  if (GFS_IS_FULL (GFS_VALUE (cell, GFS_VARIABLE (t))))
    return 0;

  guint n = 0;
  FttVector q;
  ftt_cell_pos (cell, &q);
  gdouble h = ftt_cell_size (cell);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&m->x)[c] = GFS_VALUE (cell, t->m[c]);
  gdouble alpha = GFS_VALUE (cell, t->alpha);

#if FTT_2D
  gdouble x, y;

  if (fabs (m->y) > EPS) {
    y = (alpha - m->x)/m->y;
    if (y >= 0. && y <= 1.) {
      p[n].x = q.x + h/2.; p[n].y = q.y + h*(y - 0.5); p[n++].z = 0.;
    }
  }
  if (fabs (m->x) > EPS) {
    x = (alpha - m->y)/m->x;
    if (x >= 0. && x <= 1.) {
      p[n].x = q.x + h*(x - 0.5); p[n].y = q.y + h/2.; p[n++].z = 0.;
    }
  }
  if (fabs (m->y) > EPS) {
    y = alpha/m->y;
    if (y >= 0. && y <= 1.) {
      p[n].x = q.x - h/2.; p[n].y = q.y + h*(y - 0.5); p[n++].z = 0.;
    }
  }
  if (fabs (m->x) > EPS) {
    x = alpha/m->x;
    if (x >= 0. && x <= 1.) {
      p[n].x = q.x + h*(x - 0.5); p[n].y = q.y - h/2.; p[n++].z = 0.;
    }
  }
  g_assert (n <= 2);
#else /* 3D */
  gdouble max = fabs (m->x);
  c = FTT_X;
  if (fabs (m->y) > max) {
    max = fabs (m->y);
    c = FTT_Y;
  }
  if (fabs (m->z) > max)
    c = FTT_Z;
  q.x -= h/2.; q.y -= h/2.; q.z -= h/2.;
  (&q.x)[c] += h*alpha/(&m->x)[c];
  FttVector m1 = *m;
  gts_vector_normalize (&m1.x);

  FttDirection d[12];
  n = gfs_cut_cube_vertices (cell, -1, &q, &m1, p, d, NULL, NULL);
  g_assert (n <= 6);
#endif /* 3D */
  return n;
}

/**
 * gfs_vof_facet_distance2:
 * @cell: a #FttCell.
 * @t: a #GfsVariableTracerVOF.
 * @p: a #GtsPoint.
 *
 * Returns: the square of the distance between point @p and the
 * VOF-reconstructed interface facet defined by @t or %GFS_NODATA if
 * @cell does not contain an interface.
 */
gdouble gfs_vof_facet_distance2 (FttCell * cell,
				 GfsVariableTracerVOF * t,
				 GtsPoint * p)
{
  g_return_val_if_fail (cell != NULL, GFS_NODATA);
  g_return_val_if_fail (t != NULL, GFS_NODATA);
  g_return_val_if_fail (p != NULL, GFS_NODATA);

  if (GFS_IS_FULL (GFS_VALUE (cell, GFS_VARIABLE (t))))
    return GFS_NODATA;

  FttVector q, m;
  ftt_cell_pos (cell, &q);
  gdouble h = ftt_cell_size (cell), lambda = 0., norm2 = 0.;
  FttComponent c;
  q.x -= h/2.; q.y -= h/2.; q.z -= h/2.;
  /* compute position of closest point on VOF plane m*x + m*y + m*z = alpha */
  for (c = 0; c < FTT_DIMENSION; c++) {
    (&m.x)[c] = GFS_VALUE (cell, t->m[c]);
    lambda += (&m.x)[c]*((&p->x)[c] - (&q.x)[c])/h;
    norm2 += (&m.x)[c]*(&m.x)[c];
  }
  gdouble alpha = GFS_VALUE (cell, t->alpha);
  g_assert (norm2 > 0.);
  lambda = (lambda - alpha)/norm2;

  FttVector o;
  for (c = 0; c < FTT_DIMENSION; c++) {
    (&o.x)[c] = ((&p->x)[c] - (&q.x)[c])/h - lambda*(&m.x)[c];
    if ((&o.x)[c] <= 0. || (&o.x)[c] >= 1.) {
      /* closest point on VOF plane is not within cell
	 return minimum distance from facet edges */
      FttVector q[FTT_DIMENSION*(FTT_DIMENSION - 1) + 1];
      gdouble dmin = G_MAXDOUBLE;
      guint i, n = gfs_vof_facet (cell, t, q, &m);
#if !FTT_2D
      if (n > 2)
	q[n++] = q[0];
#endif
      for (i = 0; i < n - 1; i++) {
	GtsPoint p1, p2;
	p1.x = q[i].x; p1.y = q[i].y; p1.z = q[i].z;
	p2.x = q[i + 1].x; p2.y = q[i + 1].y; p2.z = q[i + 1].z;
	GtsSegment s;
	s.v1 = (GtsVertex *) &p1; s.v2 = (GtsVertex *) &p2;
	gdouble d = gts_point_segment_distance2 (p, &s);
	if (d < dmin)
	  dmin = d;
      }
      return dmin == G_MAXDOUBLE ? GFS_NODATA : dmin;
    }
  }
  return h*h*lambda*lambda*norm2;
}

/**
 * gfs_vof_center:
 * @cell: a #FttCell.
 * @t: a #GfsVariableTracerVOF.
 * @p: a #FttVector.
 *
 * Fills @p with the coordinates of the center of mass of the
 * VOF-reconstructed interface facet defined by @t.
 *
 * Returns: the area (length in 2D) of the VOF-reconstructed facet or 0 if the
 * cell is not cut by the interface.
 */
gdouble gfs_vof_center (FttCell * cell, GfsVariableTracerVOF * t, FttVector * p)
{
  g_return_val_if_fail (cell != NULL, FALSE);
  g_return_val_if_fail (t != NULL, FALSE);
  g_return_val_if_fail (p != NULL, 0);

  if (GFS_IS_FULL (GFS_VALUE (cell, GFS_VARIABLE (t))))
    return 0.;

  FttVector m, o;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&m.x)[c] = GFS_VALUE (cell, t->m[c]);
  gdouble area = gfs_plane_area_center (&m, GFS_VALUE (cell, t->alpha), p);
  ftt_cell_pos (cell, &o);
  gdouble h = ftt_cell_size (cell);
  for (c = 0; c < FTT_DIMENSION; c++)
    (&p->x)[c] = (&o.x)[c] + h*((&p->x)[c] - 0.5);
  return area;
}

static gdouble fraction (FttVector * p,
			 guint level,
			 GfsVariable * v)
{
  FttCell * cell = gfs_domain_boundary_locate (v->domain, *p, level, NULL);
  if (cell)
    return gfs_vof_interpolate (cell, p, level, GFS_VARIABLE_TRACER_VOF (v));
  else /* fixme: boundary conditions? */
    return 2.;
}

#define NMAX 10

#define ADD_H(f) { H += f; n++; }
#define SIGN(x) ((x) > 0. ? 1. : -1.)

static gdouble local_height (FttVector * p,
			     FttVector * origin,
			     guint level,
			     GfsVariable * v,
			     FttDirection d,
			     GtsVector interface)
{
  gdouble h = ftt_level_size (level), h1 = d % 2 ? - h : h, H = 0.;
  gdouble right = fraction (p, level, v), left = right;
  FttVector pright = *p, pleft = pright;
  FttComponent c = d/2;
  guint n = 0;

  ADD_H (right);
  gboolean found_interface = (right > 0.);
  while (n < NMAX && (!found_interface || !GFS_IS_FULL (right))) {
    (&pright.x)[c] += h1;
    right = fraction (&pright, level, v);
    if (right > 1.)
      return G_MAXDOUBLE;
    ADD_H (right);
    if (!GFS_IS_FULL (right))
      found_interface = TRUE;
  }
  if (right != 1.)
    return G_MAXDOUBLE;

  found_interface = (left < 1.);
  while (n < NMAX && (!found_interface || !GFS_IS_FULL (left))) {
    (&pleft.x)[c] -= h1;
    left = fraction (&pleft, level, v);
    if (left > 1.)
      return G_MAXDOUBLE;
    ADD_H (left);
    if (!GFS_IS_FULL (left))
      found_interface = TRUE;
  }
  if (left != 0.)
    return G_MAXDOUBLE;

  H -= ((&pright.x)[c] - (&origin->x)[c])/h1 + 0.5;
  interface[0] = (p->x - origin->x)/h;
  interface[1] = (p->y - origin->y)/h;
  interface[2] = (p->z - origin->z)/h;
  interface[c] = - SIGN (h1)*H;
  return H;
}

/* fixme: does not work for periodic boundary conditions along direction c
 * for cells close to the boundary */
static gboolean curvature_along_direction (FttCell * cell, 
					   GfsVariableTracerVOF * t,
					   FttComponent c,
					   gdouble * kappa,
					   gdouble * kmax,
					   GtsVector * interface,
					   guint * n)
{
  GfsVariable * v = GFS_VARIABLE (t);

  FttVector m;
  FttComponent i;
  for (i = 0; i < FTT_DIMENSION; i++)
    (&m.x)[i] = GFS_VALUE (cell, t->m[i]);
  FttDirection d = 2*c + ((&m.x)[c] > 0.);

  FttVector p;
  ftt_cell_pos (cell, &p);
  guint level = ftt_cell_level (cell);
  gdouble size = ftt_level_size (level), H;

  gboolean found_all_heights = TRUE;
  H = local_height (&p, &p, level, v, d, interface[*n]);
  if (H == G_MAXDOUBLE)
    found_all_heights = FALSE;
  else
    (*n)++;
#if 0
  if (H < -0.5 || H > 0.5)
    found_all_heights = FALSE;
#endif

#ifdef FTT_2D
  FttComponent cp = FTT_ORTHOGONAL_COMPONENT (c);
  FttVector q = p;
  gdouble h[2];
  (&q.x)[cp] += size;
  h[0] = local_height (&q, &p, level, v, d, interface[*n]);
  if (h[0] == G_MAXDOUBLE)
    found_all_heights = FALSE;
  else
    (*n)++;

  q = p;
  (&q.x)[cp] -= size;
  h[1] = local_height (&q, &p, level, v, d, interface[*n]);
  if (h[1] == G_MAXDOUBLE)
    found_all_heights = FALSE;
  else
    (*n)++;

  if (found_all_heights) {
    gdouble hxx = h[0] - 2.*H + h[1];
    gdouble hx = (h[0] - h[1])/2.;
    gdouble dnm = 1. + hx*hx;
    *kappa = hxx/(size*sqrt (dnm*dnm*dnm));
    if (kmax)
      *kmax = fabs (*kappa);
    if (GFS_IS_AXI (v->domain)) {
      gdouble nr, r = p.y;
      if (c == FTT_X)
	nr = hx;
      else {
	r += (d == FTT_TOP ? - H*size : H*size);
	nr = (d == FTT_TOP ? 1. : -1.);
      }
      /* limit the minimum radius to half the grid size */
      gdouble kaxi = nr/MAX (sqrt(dnm)*r, size/2.);
      *kappa += kaxi;
      if (kmax)
	*kmax = MAX (*kmax, fabs (kaxi));
    }
  }
#else  /* 3D */  
  static FttComponent or[3][2] = { { FTT_Y, FTT_Z }, { FTT_X, FTT_Z }, { FTT_X, FTT_Y } };
  gdouble h[3][3];
  gint x, y;

  for (x = -1; x <= 1; x++)
    for (y = -1; y <= 1; y++)
      if (x != 0 || y != 0) {
	FttVector q = p;
	(&q.x)[or[c][0]] += size*x;
	(&q.x)[or[c][1]] += size*y;
	h[x + 1][y + 1] = local_height (&q, &p, level, v, d, interface[*n]);
	if (h[x + 1][y + 1] == G_MAXDOUBLE)
	  found_all_heights = FALSE;
	else
	  (*n)++;
      }

  if (found_all_heights) {
    gdouble hxx = h[2][1] - 2.*H + h[0][1];
    gdouble hyy = h[1][2] - 2.*H + h[1][0];
    gdouble hx = (h[2][1] - h[0][1])/2.;
    gdouble hy = (h[1][2] - h[1][0])/2.;
    gdouble hxy = (h[2][2] + h[0][0] - h[2][0] - h[0][2])/4.;
    gdouble dnm = 1. + hx*hx + hy*hy; 
    *kappa = (hxx + hyy + hxx*hy*hy + hyy*hx*hx - 2.*hxy*hx*hy)/(size*sqrt (dnm*dnm*dnm));  
    if (kmax) {
      gdouble km = *kappa/2.;
      /* Gaussian curvature */
      gdouble kg = (hxx*hyy - hxy*hxy)/(size*size*dnm*dnm);
      gdouble a = km*km - kg;
      *kmax = fabs (km);
      if (a >= 0.)
	*kmax += sqrt (a);
    }
  }
#endif /* 3D */
  return found_all_heights;
}

#define PARABOLA_FIT_CENTER_WEIGHT .1
#define PARABOLA_SIMPLER 0

typedef struct {
  GtsVector o;
#if FTT_2D /* y = a[0]*x^2 + a[0]*x + a[1] */
  GtsVector m;
  GtsMatrix * M;
  GtsVector rhs, a;
#else /* 3D */
# if PARABOLA_SIMPLER /* z = a[0]*x^2 + a[1]*y^2 + a[2]*x*y */
  GtsMatrix * M;
  GtsVector rhs, a;
# else /* z = a[0]*x^2 + a[1]*y^2 + a[2]*x*y + a[3]*x + a[4]*y + a[5] */
  gdouble ** M, rhs[6], a[6];
# endif
  GtsVector t[3];
#endif /* 3D */
} ParabolaFit;

static void parabola_fit_init (ParabolaFit * p, FttVector * o, FttVector * m)
{
  p->o[0] = o->x; p->o[1] = o->y; p->o[2] = o->z;
#if FTT_2D
  p->m[0] = m->x; p->m[1] = m->y; p->m[2] = 0.;
  gts_vector_normalize (p->m);
  p->M = gts_matrix_zero (NULL);
  p->rhs[0] = p->rhs[1] = p->rhs[2] = 0.;
#else /* 3D */
  gdouble max;
  GtsVector nx = {0., 0., 0.}, ny, nz;
  guint d = 0;

  nz[0] = m->x; nz[1] = m->y; nz[2] = m->z;
  gts_vector_normalize (nz);
  max = nz[0]*nz[0];
  /* build a vector orthogonal to nz */
  if (nz[1]*nz[1] > max) { max = nz[1]*nz[1]; d = 1; }
  if (nz[2]*nz[2] > max) d = 2;
  switch (d) {
  case 0: nx[0] = - nz[2]/nz[0]; nx[2] = 1.0; break;
  case 1: nx[1] = - nz[2]/nz[1]; nx[2] = 1.0; break;
  case 2: nx[2] = - nz[0]/nz[2]; nx[0] = 1.0; break;
  }
  gts_vector_normalize (nx);

  /* build a second vector orthogonal to nx and nz */
  gts_vector_cross (ny, nz, nx);

  /* transformation matrix from (i,j,k) to (nx, ny, nz) */
  p->t[0][0] = nx[0]; p->t[0][1] = nx[1]; p->t[0][2] = nx[2];
  p->t[1][0] = ny[0]; p->t[1][1] = ny[1]; p->t[1][2] = ny[2];
  p->t[2][0] = nz[0]; p->t[2][1] = nz[1]; p->t[2][2] = nz[2];

# if PARABOLA_SIMPLER
  p->M = gts_matrix_zero (NULL);
  p->rhs[0] = p->rhs[1] = p->rhs[2] = 0.;
# else
  p->M = gfs_matrix_new (6, 6, sizeof (gdouble));
  p->rhs[0] = p->rhs[1] = p->rhs[2] = p->rhs[3] = p->rhs[4] = p->rhs[5] = 0.;
# endif
#endif /* 3D */
}

static void parabola_fit_add (ParabolaFit * p, GtsVector m, gdouble w)
{
#if FTT_2D
  gdouble x1 = m[0] - p->o[0];
  gdouble y1 = m[1] - p->o[1];
  gdouble x = p->m[1]*x1 - p->m[0]*y1;
  gdouble y = p->m[0]*x1 + p->m[1]*y1;
  gdouble x2 = w*x*x, x3 = x2*x, x4 = x3*x;
  p->M[0][0] += x4;
  p->M[1][0] += x3; p->M[1][1] += x2;
  p->M[2][1] += w*x; p->M[2][2] += w;
  p->rhs[0] += x2*y;
  p->rhs[1] += w*x*y;
  p->rhs[2] += w*y;
#else /* 3D */
  gdouble x1 = m[0] - p->o[0];
  gdouble y1 = m[1] - p->o[1];
  gdouble z1 = m[2] - p->o[2];
  gdouble x = p->t[0][0]*x1 + p->t[0][1]*y1 + p->t[0][2]*z1;
  gdouble y = p->t[1][0]*x1 + p->t[1][1]*y1 + p->t[1][2]*z1;
  gdouble z = p->t[2][0]*x1 + p->t[2][1]*y1 + p->t[2][2]*z1;
  gdouble x2 = x*x, x3 = x2*x, x4 = x3*x;
  gdouble y2 = y*y, y3 = y2*y, y4 = y3*y;
# if PARABOLA_SIMPLER
  p->M[0][0] += w*x4;
  p->M[1][0] += w*x2*y2; p->M[1][1] += w*y4;
  p->M[2][0] += w*x3*y;  p->M[2][1] += w*x*y3;
  p->rhs[0] += w*z*x2;   p->rhs[1] += w*z*y2;  p->rhs[2] += w*z*x*y;
# else
  p->M[0][0] += w*x4; p->M[1][1] += w*y4; p->M[2][2] += w*x2*y2; 
  p->M[3][3] += w*x2; p->M[4][4] += w*y2; p->M[5][5] += w;
  p->M[0][2] += w*x3*y; p->M[0][3] += w*x3; p->M[0][4] += w*x2*y;
  p->M[1][2] += w*x*y3; p->M[1][3] += w*x*y2; p->M[1][4] += w*y3;
  p->M[2][5] += w*x*y;
  p->M[3][5] += w*x;
  p->M[4][5] += w*y;
  p->rhs[0] += w*x2*z; p->rhs[1] += w*y2*z; p->rhs[2] += w*x*y*z;
  p->rhs[3] += w*x*z; p->rhs[4] += w*y*z; p->rhs[5] += w*z;
# endif
#endif /* 3D */
}

static void parabola_fit_solve (ParabolaFit * p)
{
#if FTT_2D
  p->M[0][1] = p->M[1][0];
  p->M[0][2] = p->M[2][0] = p->M[1][1];
  p->M[1][2] = p->M[2][1];
  GtsMatrix * M = gts_matrix3_inverse ((GtsMatrix *) p->M);
  if (M) {
    p->a[0] = M[0][0]*p->rhs[0] + M[0][1]*p->rhs[1] + M[0][2]*p->rhs[2];
    p->a[1] = M[1][0]*p->rhs[0] + M[1][1]*p->rhs[1] + M[1][2]*p->rhs[2];
    gts_matrix_destroy (M);
  }
  else /* this may be a degenerate/isolated interface fragment */
    p->a[0] = p->a[1] = 0.;
#else /* 3D */
# if PARABOLA_SIMPLER
  p->M[0][1] = p->M[1][0]; p->M[0][2] = p->M[2][0];
  p->M[1][2] = p->M[2][1]; p->M[2][2] = p->M[1][0];
  GtsMatrix * M = gts_matrix3_inverse ((GtsMatrix *) p->M);
  if (M) {
    p->a[0] = M[0][0]*p->rhs[0] + M[0][1]*p->rhs[1] + M[0][2]*p->rhs[2];
    p->a[1] = M[1][0]*p->rhs[0] + M[1][1]*p->rhs[1] + M[1][2]*p->rhs[2];
    p->a[2] = M[2][0]*p->rhs[0] + M[2][1]*p->rhs[1] + M[2][2]*p->rhs[2];
    gts_matrix_destroy (M);
  }
  else /* this may be a degenerate/isolated interface fragment */
    p->a[0] = p->a[1] = p->a[2] = 0.;
# else
  p->M[0][1] = p->M[2][2]; p->M[0][5] = p->M[3][3];
  p->M[1][5] = p->M[4][4];
  p->M[2][3] = p->M[0][4]; p->M[2][4] = p->M[1][3];
  p->M[3][4] = p->M[2][5];
  guint i, j;
  for (i = 1; i < 6; i++)
    for (j = 0; j < i; j++)
      p->M[i][j] = p->M[j][i];
  if (gfs_matrix_inverse (p->M, 6, 1e-10)) {
    for (i = 0; i < 6; i++) {
      p->a[i] = 0.;
      for (j = 0; j < 6; j++)
	p->a[i] += p->M[i][j]*p->rhs[j];
    }
  }
  else { /* this may be a degenerate/isolated interface fragment */
    g_warning ("singular matrix");
    p->a[0] = p->a[1] = p->a[2] = 0.;
  }
# endif
#endif /* 3D */
}

static gdouble parabola_fit_curvature (ParabolaFit * p, gdouble kappamax,
				       gdouble * kmax)
{
  gdouble kappa;
#if FTT_2D
  gdouble dnm = 1. + p->a[1]*p->a[1];
  kappa = 2.*p->a[0]/sqrt (dnm*dnm*dnm);
  if (kmax)
    *kmax = fabs (kappa);
#else /* 3D */
  gdouble hxx = 2.*p->a[0];
  gdouble hyy = 2.*p->a[1];
  gdouble hxy = p->a[2];
  gdouble hx, hy;
# if PARABOLA_SIMPLER
  hx = hy = 0.;
# else
  hx = p->a[3];
  hy = p->a[4];
# endif
  gdouble dnm = 1. + hx*hx + hy*hy;
  kappa = (hxx + hyy + hxx*hy*hy + hyy*hx*hx - 2.*hxy*hx*hy)/sqrt (dnm*dnm*dnm);
  if (kmax) {
    gdouble kg = (hxx*hyy - hxy*hxy)/(dnm*dnm);
    gdouble a = kappa*kappa/4. - kg;
    *kmax = fabs (kappa/2.);
    if (a >= 0.)
      *kmax += sqrt (a);
  }
#endif /* 3D */
  if (fabs (kappa) > kappamax) {
    if (kmax)
      *kmax = kappamax;
    return kappa > 0. ? kappamax : - kappamax;
  }
  return kappa;
}

#if FTT_2D
static void parabola_fit_axi_curvature (const ParabolaFit * p, gdouble r, gdouble h,
					gdouble * kappa, gdouble * kmax)
{
  gdouble nr = (p->m[0]*p->a[1] + p->m[1])/sqrt (1. + p->a[1]*p->a[1]);
  gdouble kaxi = - nr/MAX(r, h/2.); /* limit the minimum radius to half the grid size */
  *kappa += kaxi;
  if (kmax)
    *kmax = MAX (*kmax, fabs (kaxi));
}
#endif /* 2D */

static void parabola_fit_destroy (ParabolaFit * p)
{
#if (FTT_2D || PARABOLA_SIMPLER)
  gts_matrix_destroy (p->M);
#else
  gfs_matrix_free (p->M);
#endif
}

static void add_vof_center (FttCell * cell, FttVector * p, guint level,
			    FttVector * origin,
			    GfsVariableTracerVOF * t,
			    ParabolaFit * fit, gdouble w)
{
  gdouble f = GFS_VALUE (cell, GFS_VARIABLE (t));
  if (!GFS_IS_FULL (f)) {
    FttVector m, c;
    gdouble alpha = gfs_vof_plane_interpolate (cell, p, level, t, &m);
    gdouble area = gfs_plane_area_center (&m, alpha, &c);
    gdouble h = ftt_level_size (level);
    FttComponent i;
    for (i = 0; i < FTT_DIMENSION; i++)
      (&c.x)[i] = ((&p->x)[i] - (&origin->x)[i])/h + (&c.x)[i] - 0.5;
    parabola_fit_add (fit, &c.x, w*area);
  }
}

static void fit_from_fractions (FttCell * cell, GfsVariable * v, ParabolaFit * fit)
{
  gdouble h = ftt_cell_size (cell);
  guint level = ftt_cell_level (cell);
  gint x, y, z = 0;
  FttVector p;
  
  ftt_cell_pos (cell, &p);
#if !FTT_2D
  for (z = -1; z <= 1; z++)
#endif
    for (x = -1; x <= 1; x++)
      for (y = -1; y <= 1; y++)
	if (x != 0 || y != 0 || z != 0) {
	  FttVector o;
	  o.x = p.x + h*x; o.y = p.y + h*y; o.z = p.z + h*z;
	  FttCell * neighbor = gfs_domain_boundary_locate (v->domain, o, level, NULL);
	  if (neighbor)
	    add_vof_center (neighbor, &o, level, &p, GFS_VARIABLE_TRACER_VOF (v),
			    fit, 1.);
	}
}

/**
 * gfs_fit_curvature:
 * @cell: a #FttCell containing an interface.
 * @t: a #GfsVariableTracerVOF.
 * @kmax: a pointer or %NULL. 
 *
 * Computes an approximation of the curvature of the interface
 * contained in @cell using paraboloid fitting of the centroids of the
 * reconstructed interface segments.
 *
 * If @kmax is not %NULL, it is filled with the absolute value of the
 * maximum surface curvature (note that in 2D this is just the absolute value of
 * the mean curvature).
 *
 * Returns: (double in 3D) the mean curvature of the interface contained in @cell.
 */
gdouble gfs_fit_curvature (FttCell * cell, GfsVariableTracerVOF * t, gdouble * kmax)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (t != NULL, 0.);

  GfsVariable * v = GFS_VARIABLE (t);
  g_return_val_if_fail (!GFS_IS_FULL (GFS_VALUE (cell,  v)), 0.);

  FttVector m;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&m.x)[c] = GFS_VALUE (cell, t->m[c]);

  ParabolaFit fit;
  FttVector p, fc;
  ftt_cell_pos (cell, &p);
  gdouble area = gfs_vof_center (cell, t, &fc);
  gdouble h = ftt_cell_size (cell);
  fc.x = (fc.x - p.x)/h;
  fc.y = (fc.y - p.y)/h;
  fc.z = (fc.z - p.z)/h;
  parabola_fit_init (&fit, &fc, &m);
  parabola_fit_add (&fit, &fc.x, area);
  fit_from_fractions (cell, GFS_VARIABLE (t), &fit);
  parabola_fit_solve (&fit);
  gdouble kappa = parabola_fit_curvature (&fit, 2., kmax)/h;
  if (kmax)
    *kmax /= h;
#if FTT_2D
  if (GFS_IS_AXI (v->domain))
    parabola_fit_axi_curvature (&fit, fc.y*h + p.y, h, &kappa, kmax);
#endif
  parabola_fit_destroy (&fit);
  return kappa;
}

#if FTT_2D
# define NI 3
#else
# define NI 9
#endif

static void orientation (FttVector * m, FttComponent * c)
{
  FttComponent i, j;
  for (i = 0; i < FTT_DIMENSION; i++)
    c[i] = i;
  for (i = 0; i < FTT_DIMENSION - 1; i++)
    for (j = 0; j < FTT_DIMENSION - 1 - i; j++)
      if (fabs ((&m->x)[c[j + 1]]) > fabs ((&m->x)[c[j]])) {
	FttComponent tmp = c[j];
	c[j] = c[j + 1];
	c[j + 1] = tmp;
      }
}

static guint independent_positions (GtsVector * interface, guint n)
{
  if (n < 2)
    return n;

  guint j, ni = 1;
  for (j = 1; j < n; j++) {
    guint i;
    gboolean depends = FALSE;
    for (i = 0; i < j && !depends; i++) {
      gdouble d2 = 0.;
      FttComponent c;
      for (c = 0; c < FTT_DIMENSION; c++)
	d2 += (interface[i][c] - interface[j][c])*(interface[i][c] - interface[j][c]);
      depends = d2 < 0.5*0.5;
    }
    ni += !depends;
  }
  return ni;
}

/**
 * gfs_height_curvature:
 * @cell: a #FttCell containing an interface.
 * @t: a #GfsVariableTracerVOF.
 * @kmax: a pointer or %NULL.
 *
 * An implementation of the Height-Function (HF) method generalised to
 * adaptive meshes.
 *
 * If @kmax is not %NULL, it is filled with the absolute value of the
 * maximum surface curvature (note that in 2D this is just the
 * absolute value of the mean curvature).
 *
 * Returns: (double in 3D) the mean curvature of the interface
 * contained in @cell, or %GFS_NODATA if the HF method could not
 * compute a consistent curvature.
 */
gdouble gfs_height_curvature (FttCell * cell, GfsVariableTracerVOF * t, gdouble * kmax)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (t != NULL, 0.);

  GfsVariable * v = GFS_VARIABLE (t);
  gdouble f = GFS_VALUE (cell,  v);
  g_return_val_if_fail (!GFS_IS_FULL (f), 0.);

  FttVector m;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&m.x)[c] = GFS_VALUE (cell, t->m[c]);

  FttComponent try[FTT_DIMENSION];
  orientation (&m, try); /* sort directions according to normal */

  gdouble kappa = 0.;
  GtsVector interface[FTT_DIMENSION*NI];
  guint n = 0;
  for (c = 0; c < FTT_DIMENSION; c++) /* try each direction */
    if (curvature_along_direction (cell, t, try[c], &kappa, kmax, interface, &n))
      return kappa;

  /* Could not compute curvature from the simple algorithm along any direction:
   * Try parabola fitting of the collected interface positions */

  if (independent_positions (interface, n) < 3*(FTT_DIMENSION - 1)) {
    if (kmax)
      *kmax = GFS_NODATA;
    return GFS_NODATA;
  }

  gdouble h = ftt_cell_size (cell);
  ParabolaFit fit;
  guint j;
  
  FttVector p, fc;
  ftt_cell_pos (cell, &p);
  gdouble area = gfs_vof_center (cell, t, &fc);
  fc.x = (fc.x - p.x)/h;
  fc.y = (fc.y - p.y)/h;
  fc.z = (fc.z - p.z)/h;
  parabola_fit_init (&fit, &fc, &m);
#if FTT_2D
  parabola_fit_add (&fit, &fc.x, PARABOLA_FIT_CENTER_WEIGHT);
#elif !PARABOLA_SIMPLER
  parabola_fit_add (&fit, &fc.x, area*100.);
#endif
  for (j = 0; j < n; j++)
    parabola_fit_add (&fit, interface[j], 1.);
  parabola_fit_solve (&fit);
  kappa = parabola_fit_curvature (&fit, 2., kmax)/h;
  if (kmax)
    *kmax /= h;
#if FTT_2D
  if (GFS_IS_AXI (v->domain))
    parabola_fit_axi_curvature (&fit, fc.y*h + p.y, h, &kappa, kmax);
#endif
  parabola_fit_destroy (&fit);
  return kappa;
}

/* Returns: the height @h of the neighboring column in direction @d or
   GFS_NODATA if it is undefined. Also fills @x with the coordinates
   of the cell (3/4, 1 or 3/2 depending on its relative level). */
static gdouble neighboring_column (FttCell * cell, 
				   GfsVariable * h, FttComponent c, gdouble orientation,
				   FttDirection d, gdouble * x)
{
  FttCell * n = ftt_cell_neighbor (cell, d);
  if (!n)
    return GFS_NODATA;
  if (ftt_cell_level (cell) == ftt_cell_level (n)) {
    if (GFS_HAS_DATA (n, h)) {
      *x = 1.;
      return GFS_VALUE (n, h);
    }
    else if (FTT_CELL_IS_LEAF (n))
      return GFS_NODATA;
    /* check finer neighbors */
    FttCellChildren child;
    int i, m = ftt_cell_children_direction (n, FTT_OPPOSITE_DIRECTION (d), &child);
    for (i = 0; i < m; i++)
      if (child.c[i] && GFS_HAS_DATA (child.c[i], h)) {
	FttVector p;
	ftt_cell_relative_pos (child.c[i], &p);
	*x = 3./4.;
	return GFS_VALUE (child.c[i], h)/2. + orientation*(&p.x)[c];
      }
    return GFS_NODATA;
  }
  else if (GFS_HAS_DATA (n, h)) {
    /* coarser neighbor */
    FttVector p;
    ftt_cell_relative_pos (cell, &p);
    *x = 3./2.;
    return (GFS_VALUE (n, h) - orientation*(&p.x)[c])*2.;
  }
  return GFS_NODATA;
}

static void curvature_from_h (FttCell * cell, GfsDomain * domain,
			      gdouble x[3], gdouble h[3],
			      gdouble orientation, FttComponent c,
			      gdouble * kappa, gdouble * kmax)
{
  gdouble size = ftt_cell_size (cell);
  gdouble det = x[0]*x[1]*(x[0] - x[1]), a = x[1]*(h[0] - h[2]), b = x[0]*(h[1] - h[2]);
  gdouble hxx = 2.*(a - b)/det;
  gdouble hx = (x[0]*b - x[1]*a)/det;
  gdouble dnm = 1. + hx*hx;
  *kappa = hxx/(size*sqrt (dnm*dnm*dnm));
  if (kmax)
    *kmax = fabs (*kappa);
  if (GFS_IS_AXI (domain)) {
    FttVector p;
    ftt_cell_pos (cell, &p);
    gdouble nr, r = p.y;
    if (c == FTT_X)
      nr = hx;
    else {
      r += orientation*h[2]*size;
      nr = - orientation;
    }
    /* limit the minimum radius to half the grid size */
    gdouble kaxi = nr/MAX (sqrt(dnm)*r, size/2.);
    *kappa += kaxi;
    if (kmax)
      *kmax = MAX (*kmax, fabs (kaxi));
  }
}

/**
 * gfs_curvature_along_direction:
 * @cell: a #FttCell.
 * @t: a #GfsVariableTracerVOFHeight.
 * @c: x, y or z.
 * @kappa: the curvature.
 * @kmax: the maximum curvature.
 *
 * Tries to compute an interface curvature for @cell using
 * height-functions on equally-spaced columns in direction @c.
 *
 * Returns: %TRUE if the curvature was successfully computed, %FALSE
 * otherwise.
 */
gboolean gfs_curvature_along_direction (FttCell * cell, 
					GfsVariableTracerVOFHeight * t,
					FttComponent c,
					gdouble * kappa,
					gdouble * kmax)
{
  g_return_val_if_fail (cell != NULL, FALSE);
  g_return_val_if_fail (t != NULL, FALSE);
  g_return_val_if_fail (kappa != NULL, FALSE);

#ifdef FTT_2D
  gdouble orientation;
  GfsVariable * hv = gfs_closest_height (cell, t, c, &orientation);
  if (!hv)
    return FALSE;
  else if (fabs (GFS_VALUE (cell, hv)) > 1.)
    return FALSE; /* interface is too far */

  FttComponent oc = FTT_ORTHOGONAL_COMPONENT (c);
  gdouble x[3], h[3];
  h[2] = GFS_VALUE (cell, hv); x[2] = 0.;
  h[0] = neighboring_column (cell, hv, c, orientation, 2*oc, &x[0]);
  if (h[0] != GFS_NODATA && x[0] == 1.) {
    h[1] = neighboring_column (cell, hv, c, orientation, 2*oc + 1, &x[1]);
    if (h[1] != GFS_NODATA && x[1] == 1.) {
      x[1] = - x[1];
      curvature_from_h (cell, GFS_VARIABLE (t)->domain, x, h, orientation, c, kappa, kmax);
      return TRUE;
    }
  }
#else /* 3D */
  g_assert_not_implemented ();
#endif /* 3D */

  return FALSE;
}

static gboolean curvature_along_direction_new (FttCell * cell, 
					       GfsVariableTracerVOFHeight * t,
					       FttComponent c,
					       gdouble * kappa,
					       gdouble * kmax,
					       GtsVector * interface,
					       guint * nb)
{
#ifdef FTT_2D
  gdouble orientation;
  GfsVariable * hv = gfs_closest_height (cell, t, c, &orientation);
  FttComponent oc = FTT_ORTHOGONAL_COMPONENT (c);
  if (!hv) {
    /* no data for either directions, look "right" and "left" to
       collect potential interface positions */
    hv = gfs_closest_height (ftt_cell_neighbor (cell, 2*oc), t, c, &orientation);
    if (!hv)
      hv = gfs_closest_height (ftt_cell_neighbor (cell, 2*oc + 1), t, c, &orientation);
    if (!hv) /* give up */
      return FALSE;
  }
  else if (fabs (GFS_VALUE (cell, hv)) > 1.)
    return FALSE; /* interface is too far */

  gdouble x[3], h[3];
  h[2] = GFS_VALUE (cell, hv); x[2] = 0.;
  h[0] = neighboring_column (cell, hv, c, orientation, 2*oc, &x[0]);
  h[1] = neighboring_column (cell, hv, c, orientation, 2*oc + 1, &x[1]);
  x[1] = - x[1];

  if (h[2] != GFS_NODATA && h[0] != GFS_NODATA && h[1] != GFS_NODATA) {
    curvature_from_h (cell, GFS_VARIABLE (t)->domain, x, h, orientation, c, kappa, kmax);
    return TRUE;
  }
  else { /* h[2] == GFS_NODATA || h[0] == GFS_NODATA || h[1] == GFS_NODATA */
    /* collect interface positions (based on height function) */
    int i;
    for (i = 0; i < 3; i++)
      if (h[i] != GFS_NODATA) {
	interface[*nb][oc] = x[i];
	interface[(*nb)++][c] = orientation*h[i]; 
      }
    return FALSE;
  }
#else /* 3D */
  g_assert_not_implemented ();
#endif /* 3D */

  return FALSE;
}

/**
 * gfs_height_curvature_new:
 * @cell: a #FttCell containing an interface.
 * @t: a #GfsVariableTracerVOFHeight.
 * @kmax: a pointer or %NULL.
 *
 * Tries to estimate the curvature of an interface using
 * height-functions, either on equally-spaced columns, non-equally
 * spaced columns or using parabola fits of interface positions
 * defined using the height-functions in all directions.
 *
 * If @kmax is not %NULL, it is filled with the absolute value of the
 * maximum surface curvature (note that in 2D this is just the
 * absolute value of the mean curvature).
 *
 * Returns: (double in 3D) the mean curvature of the interface
 * contained in @cell, or %GFS_NODATA if the HF method could not
 * compute a consistent curvature.
 */
gdouble gfs_height_curvature_new (FttCell * cell, GfsVariableTracerVOFHeight * t, 
				  gdouble * kmax)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (t != NULL, 0.);

  GfsVariable * v = GFS_VARIABLE (t);
  gdouble f = GFS_VALUE (cell,  v);
  g_return_val_if_fail (!GFS_IS_FULL (f), 0.);

  FttVector m;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&m.x)[c] = GFS_VALUE (cell, GFS_VARIABLE_TRACER_VOF (t)->m[c]);

  FttComponent try[FTT_DIMENSION];
  orientation (&m, try); /* sort directions according to normal */

  gdouble kappa = 0.;
  GtsVector interface[FTT_DIMENSION*NI];
  guint n = 0;
  for (c = 0; c < FTT_DIMENSION; c++) /* try each direction */
    if (curvature_along_direction_new (cell, t, try[c], &kappa, kmax, interface, &n))
      return kappa;

  /* Could not compute curvature from the simple algorithm along any direction:
   * Try parabola fitting of the collected interface positions */

  if (independent_positions (interface, n) < 3*(FTT_DIMENSION - 1))
    return GFS_NODATA;

  gdouble h = ftt_cell_size (cell);
  ParabolaFit fit;
  guint j;
  
  FttVector p, fc;
  ftt_cell_pos (cell, &p);
  gdouble area = gfs_vof_center (cell, GFS_VARIABLE_TRACER_VOF (t), &fc);
  fc.x = (fc.x - p.x)/h;
  fc.y = (fc.y - p.y)/h;
  fc.z = (fc.z - p.z)/h;
  parabola_fit_init (&fit, &fc, &m);
#if FTT_2D
  parabola_fit_add (&fit, &fc.x, PARABOLA_FIT_CENTER_WEIGHT);
#elif !PARABOLA_SIMPLER
  parabola_fit_add (&fit, &fc.x, area*100.);
#endif
  for (j = 0; j < n; j++)
    parabola_fit_add (&fit, interface[j], 1.);
  parabola_fit_solve (&fit);
  kappa = parabola_fit_curvature (&fit, 2., kmax)/h;
  if (kmax)
    *kmax /= h;
#if FTT_2D
  if (GFS_IS_AXI (v->domain))
    parabola_fit_axi_curvature (&fit, fc.y*h + p.y, h, &kappa, kmax);
#endif
  parabola_fit_destroy (&fit);
  return kappa;
}

/**
 * gfs_vof_correctness:
 * @cell: a #FttCell.
 * @t: a #GfsVariableTracerVOF.
 *
 * An implementation of the criterion of Cerne, Petelin, Tiselj
 * (2002), to measure how well an interface is represented by a local
 * VOF field.
 *
 * Returns: the "correctness" of the interface representation.
 */
gdouble gfs_vof_correctness (FttCell * cell, GfsVariableTracerVOF * t)
{
  GfsVariable * v = GFS_VARIABLE (t);
  gdouble F(3,3,3);
  
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (t != NULL, 0.);

  if (GFS_VALUE (cell, v) <= 0. || GFS_VALUE (cell, v) >= 1.)
    return 1.;

  stencil (cell, v, f);
#if FTT_2D
  gdouble dx = f[2][0] + f[2][1] + f[2][2] - f[0][0] - f[0][1] - f[0][2];
  gdouble dy = f[0][2] + f[1][2] + f[2][2] - f[0][0] - f[1][0] - f[2][0];
  return sqrt ((dx*dx + dy*dy)/9.);
#else
  gdouble dx = (f[2][0][0] + f[2][1][0] + f[2][2][0] - f[0][0][0] - f[0][1][0] - f[0][2][0] +
		f[2][0][1] + f[2][1][1] + f[2][2][1] - f[0][0][1] - f[0][1][1] - f[0][2][1] +
		f[2][0][2] + f[2][1][2] + f[2][2][2] - f[0][0][2] - f[0][1][2] - f[0][2][2]);
  gdouble dy = (f[0][2][0] + f[1][2][0] + f[2][2][0] - f[0][0][0] - f[1][0][0] - f[2][0][0] +
		f[0][2][1] + f[1][2][1] + f[2][2][1] - f[0][0][1] - f[1][0][1] - f[2][0][1] +
		f[0][2][2] + f[1][2][2] + f[2][2][2] - f[0][0][2] - f[1][0][2] - f[2][0][2]);
  gdouble dz = (f[0][0][2] + f[1][0][2] + f[2][0][2] - f[0][0][0] - f[1][0][0] - f[2][0][0] +
		f[0][1][2] + f[1][1][2] + f[2][1][2] - f[0][1][0] - f[1][1][0] - f[2][1][0] +
		f[0][2][2] + f[1][2][2] + f[2][2][2] - f[0][2][0] - f[1][2][0] - f[2][2][0]);  
  return sqrt ((dx*dx + dy*dy + dz*dz)/27.);
#endif
}

/** \endobject{GfsVariableTracerVOF} */

/** \beginobject{GfsVariableVOFConcentration} */

static void variable_vof_concentration_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_vof_concentration_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (name)");
    return;
  }
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  GfsVariable * v = gfs_variable_from_name (domain->variables, fp->token->str);  
  if (v == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  if (!GFS_IS_VARIABLE_TRACER_VOF (v)) {
    gts_file_error (fp, "variable `%s' is not a VOF tracer", fp->token->str);
    return;
  }
  GfsVariableTracerVOF * vof = GFS_VARIABLE_TRACER_VOF (v);
  GFS_VARIABLE_VOF_CONCENTRATION (*o)->vof = vof;
  gts_container_add (GTS_CONTAINER (vof->concentrations), GTS_CONTAINEE (*o));
  gts_file_next_token (fp);
}

static void variable_vof_concentration_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_vof_concentration_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s", GFS_VARIABLE (GFS_VARIABLE_VOF_CONCENTRATION (o)->vof)->name);
}

static void variable_vof_concentration_class_init (GtsObjectClass * klass)
{
  klass->read = variable_vof_concentration_read;
  klass->write = variable_vof_concentration_write;
}

static void variable_vof_concentration_init (GfsVariable * v)
{
  /* this is taken care of by the associated VOF tracer */
  v->fine_coarse = v->coarse_fine = no_coarse_fine;
  GFS_VARIABLE_TRACER (v)->advection.scheme = GFS_NONE;
}

GfsVariableClass * gfs_variable_vof_concentration_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsVariableVOFConcentration",
      sizeof (GfsVariableVOFConcentration),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_vof_concentration_class_init,
      (GtsObjectInitFunc) variable_vof_concentration_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_tracer_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsVariableVOFConcentration} */

/** \beginobject{GfsVariableTracerVOFHeight} */

typedef struct {
  GfsVariable * f; /* volume fraction */
  GfsVariable * hb, * ht; /* heights in either orientation */
  GfsBc * angle; /* contact angle BC */
  FttComponent c; /* x, y or z */
  FttDirection d;
} HFState;

static gboolean is_interfacial (FttCell * cell, gpointer data)
{
  GfsVariable * f = data;
  return !GFS_IS_FULL (GFS_VALUE (cell, f));
}

static void undefined_height (FttCell * cell, HFState * hf)
{
  GFS_VALUE (cell, hf->hb) = GFS_NODATA;
  GFS_VALUE (cell, hf->ht) = GFS_NODATA;
}

#define HMAX 5
#define SIGN(x) ((x) > 0. ? 1. : -1.)
#define BOUNDARY_HIT (2.*HMAX)

static gint children_half_height (FttCell * cell, FttDirection d, GfsVariable * fv,
				  gdouble * H, gint * n)
{
  FttCellChildren child;
  guint i, m = ftt_cell_children_direction (cell, FTT_OPPOSITE_DIRECTION (d), &child);
  gint s = 0;
  for (i = 0; i < m; i++)
    if (child.c[i]) {
      gdouble f = GFS_VALUE (child.c[i], fv);
      if (f > 0. && f < 1.) {
	s = 0;
	break;
      }
      s = SIGN (f - 0.5);
    }
  if (s != 0)
    return s;

  *H += GFS_VALUE (cell, fv);
  (*n)++;

  m = ftt_cell_children_direction (cell, d, &child);
  for (i = 0; i < m; i++)
    if (child.c[i]) {
      gdouble f = GFS_VALUE (child.c[i], fv);
      if (f > 0. && f < 1.) {
	return 0;
      }
      s = SIGN (f - 0.5);
    }
  return s;
}

static gint half_height (FttCell * cell, GfsVariable * fv, FttDirection d,
			 gdouble * H, gint * n)
{
  gint s = 0;
  *n = 0;
  guint level = ftt_cell_level (cell);
  FttCell * neighbor = ftt_cell_neighbor (cell, d);
  while (*n < HMAX && !s && neighbor) {
    gdouble f = GFS_VALUE (neighbor, fv);
    if (f > 0. && f < 1.) { /* interfacial cell */
      if (ftt_cell_level (neighbor) < level) /* neighbor is coarser */
	return 0;
      if (GFS_CELL_IS_BOUNDARY (neighbor))
	return 2;
      if (FTT_CELL_IS_LEAF (neighbor)) {
	*H += f;
	(*n)++;
      }
      else
	s = children_half_height (neighbor, d, fv, H, n);
    }
    else /* full or empty cell */
      s = SIGN (f - 0.5);
    neighbor = ftt_cell_neighbor (neighbor, d);
  }
  return s;
}

#define DMAX 2.

static void height_propagation (FttCell * cell, HFState * hf, GfsVariable * h, gdouble orientation)
{
  guint level = ftt_cell_level (cell);
  FttDirection d;
  for (d = 2*hf->c; d <= 2*hf->c + 1; d++, orientation = - orientation) {
    gdouble H = GFS_VALUE (cell, h);
    FttCell * neighbor = ftt_cell_neighbor (cell, d);
    gboolean interface = FALSE;
    while (fabs (H) < DMAX - 1. && neighbor && !interface && 
	   ftt_cell_level (neighbor) == level) {
      H -= orientation;
      GFS_VALUE (neighbor, h) = H;
      interface = is_interfacial (neighbor, hf->f);
      neighbor = ftt_cell_neighbor (neighbor, d);
    }
  }
}

static void height (FttCell * cell, HFState * hf)
{
  if (!FTT_CELL_IS_LEAF (cell)) {
    /* look for childrens' HF */
    FttComponent c = FTT_ORTHOGONAL_COMPONENT (hf->c);
    gdouble H = 0., orientation = 0.;
    guint nc = 0;
    GfsVariable * h = NULL;
    FttDirection d;
    for (d = 2*c; d <= 2*c + 1; d++) {
      FttCellChildren child;
      int i, n = ftt_cell_children_direction (cell, d, &child);
      for (i = 0; i < n; i++)
	if (child.c[i]) {
	  if (h == NULL)
	    h = gfs_closest_height (child.c[i], 
				    GFS_VARIABLE_TRACER_VOF_HEIGHT (hf->f), hf->c, &orientation);
	  if (h != NULL && GFS_HAS_DATA (child.c[i], h)) {
	    FttVector p;
	    ftt_cell_relative_pos (child.c[i], &p);
	    H += GFS_VALUE (child.c[i], h)/2. + orientation*(&p.x)[hf->c];
	    nc++;
	    break;
	  }
	}
    }
    if (nc == 2) {
      GFS_VALUE (cell, h) = H/nc;
      height_propagation (cell, hf, h, orientation);
      return;
    }
    /* try the standard way just in case */
  }

  gdouble H = GFS_VALUE (cell, hf->f);

  /* top part of the column */
  gint nt, st = half_height (cell, hf->f, 2*hf->c, &H, &nt);
  if (!st) /* still an interfacial cell (or coarser neighboring cell found) */
    return;

  /* bottom part of the column */
  gint nb, sb = half_height (cell, hf->f, 2*hf->c + 1, &H, &nb);
  if (!sb) /* still an interfacial cell (or coarser neighboring cell found) */
    return;

  if (sb != 2 && st != 2) {
    if (st*sb > 0) /* the column does not cross the interface */
      return;
  }
  else { /* column hit a boundary */
    if (sb == 2 && st == 2) /* cannot hit a boundary on both sides */
      return;
    if (sb == 2)
      sb = - SIGN (st);
    H += BOUNDARY_HIT;
  }

  if (sb > 0) {
    GFS_VALUE (cell, hf->hb) = H - 0.5 - nb;
    height_propagation (cell, hf, hf->hb, 1.);
  }
  else {
    GFS_VALUE (cell, hf->ht) = H - 0.5 - nt;
    height_propagation (cell, hf, hf->ht, -1.);
  }
}

static GfsVariable * boundary_hit (FttCell * cell, HFState * hf)
{
  if (GFS_HAS_DATA (cell, hf->hb) && GFS_VALUE (cell, hf->hb) > BOUNDARY_HIT/2.)
    return hf->hb;
  if (GFS_HAS_DATA (cell, hf->ht) && GFS_VALUE (cell, hf->ht) > BOUNDARY_HIT/2.)
    return hf->ht;
  return NULL;
}

static void height_propagation_from_boundary (FttCell * cell, HFState * hf, GfsVariable * h)
{
  guint level = ftt_cell_level (cell);
  FttDirection d = FTT_OPPOSITE_DIRECTION (hf->d);
  gdouble orientation = (hf->d % 2 ? -1 : 1)*(h == hf->hb ? 1 : -1);
  gdouble H = GFS_VALUE (cell, h);
  cell = ftt_cell_neighbor (cell, d);
  while (cell && GFS_HAS_DATA (cell, h) && GFS_VALUE (cell, h) > BOUNDARY_HIT/2. &&
	 ftt_cell_level (cell) == level) {
    H += orientation;
    GFS_VALUE (cell, h) = H;
    cell = ftt_cell_neighbor (cell, d);
  }
  /* propagate to non-interfacial cells up to DMAX */
  while (fabs (H) < DMAX - 1. && cell && !is_interfacial (cell, hf->f) && 
	 ftt_cell_level (cell) == level) {
    H += orientation;
    GFS_VALUE (cell, h) = H;
    cell = ftt_cell_neighbor (cell, d);
  }
}

static void height_periodic_bc (FttCell * cell, HFState * hf)
{
  FttCell * neighbor = ftt_cell_neighbor (cell, hf->d);
  g_assert (GFS_CELL_IS_BOUNDARY (neighbor));
  GfsVariable * h = boundary_hit (cell, hf);
  if (h) {
    /* column hit boundary */
    GfsVariable * hn = boundary_hit (neighbor, hf);
    if (hn == h) {
      /* the column crosses the interface */
      /* propagate column height correction from one side (or PE) to the other */
      gdouble orientation = (hf->d % 2 ? -1 : 1)*(h == hf->hb ? 1 : -1);
      gdouble Hn = GFS_VALUE (neighbor, h) + 0.5 + 
	(orientation - 1.)/2. -
	2.*BOUNDARY_HIT;
      GFS_VALUE (cell, h) += Hn;
      height_propagation_from_boundary (cell, hf, h);
    }
    else {
      /* the column does not cross the interface */
      guint level = ftt_cell_level (cell);
      while (cell && GFS_HAS_DATA (cell, h) && GFS_VALUE (cell, h) > BOUNDARY_HIT/2. &&
	     ftt_cell_level (cell) == level) {
	undefined_height (cell, hf);
	cell = ftt_cell_neighbor (cell, FTT_OPPOSITE_DIRECTION (hf->d));
      }
    }
  }
  else {
    /* column did not hit a boundary, propagate height across PE boundary */
    if (GFS_HAS_DATA (neighbor, hf->hb))
      height_propagation (neighbor, hf, hf->hb, 1.);
    if (GFS_HAS_DATA (neighbor, hf->ht))
      height_propagation (neighbor, hf, hf->ht, -1.);
  }
}

#define SLOPE_MAX (2.*HMAX/3.)
#define THETA_MIN (atan(1./SLOPE_MAX))

static gdouble contact_angle_bc (FttCell * cell, HFState * hf)
{
  if (hf->angle) {
    FttCellFace f = ftt_cell_face (cell, hf->d);
    g_assert (GFS_CELL_IS_BOUNDARY (f.neighbor));
    return gfs_function_face_value (GFS_BC_VALUE (hf->angle)->val, &f)*M_PI/180.;
  }
  else
    return M_PI/2.;
}

static void height_contact_normal_bc (FttCell * cell, HFState * hf)
{
  GfsVariable * h = boundary_hit (cell, hf);
  if (h) {
    /* column hit boundary */
    FttComponent oc = FTT_ORTHOGONAL_COMPONENT (hf->c);
    FttDirection nd = 2*oc;
    FttCell * n1 = ftt_cell_neighbor (cell, nd), * n2 = ftt_cell_neighbor (cell, nd + 1);
    if (n1 && GFS_HAS_DATA (n1, h) && n2 && GFS_HAS_DATA (n2, h)) {
      /* columns are defined on either side of @cell => it is not a contact line */
      /* the boundary is dry or wet i.e. the column height is well defined */
      /* reset the BOUNDARY_HIT flag on the whole column */
      GFS_VALUE (cell, h) -= BOUNDARY_HIT;
      height_propagation_from_boundary (cell, hf, h);
    }
    else {
      GfsVariable * hb, * ht; /* use symmetries */
      if (hf->d % 2) {
	hb = hf->hb; ht = hf->ht;
      }
      else {
	hb = hf->ht; ht = hf->hb;
      }
      gdouble full_or_empty = (h != hb);
      if (n1 && GFS_VALUE (n1, hf->f) != full_or_empty) {
	n1 = n2; nd++;
      }
      if (!n1 || GFS_VALUE (n1, hf->f) != full_or_empty) {
	/* column is not neighbouring a full or empty cell => it is not a contact line */
	/* the boundary is dry or wet i.e. the column height is well defined */
	/* reset the BOUNDARY_HIT flag on the whole column */
	GFS_VALUE (cell, h) -= BOUNDARY_HIT;
	height_propagation_from_boundary (cell, hf, h);
      }
      /* contact line */
      else {
	gdouble theta = contact_angle_bc (cell, hf);
	if ((h == hb && theta < atan (SLOPE_MAX)) || 
	    (h == ht && theta > M_PI - atan (SLOPE_MAX))) {
	  gdouble orientation = (h == hb ? 1. : -1.);
	  FttVector m = { orientation*sin(theta), cos(theta), 0. };
	  gdouble alpha = gfs_plane_alpha (&m, GFS_VALUE (cell, hf->f));
	  GFS_VALUE (cell, h) = orientation*((alpha - m.x/2.)/m.y - 0.5);
	  height_propagation_from_boundary (cell, hf, h);
	  /* set height of neighbouring (non-interfacial) column */
	  /* the line below ensures that the interface does not enter
	     the non-interfacial neighbour */
	  if (orientation*alpha > orientation*m.x) alpha = m.x;
	  GFS_VALUE (n1, h) = ftt_cell_level (n1) == ftt_cell_level (cell) ? 
	    orientation*((alpha - m.x*3./2.)/m.y - 0.5) : /* neighbour at same level */
	    orientation*((alpha - m.x*2.)/m.y - 1.)/2.;   /* coarser neighbor */
	  height_propagation_from_boundary (n1, hf, h);
	}
      }
    }
  }
}

static void contact_angle_height (FttCell * cell, GfsVariable * h, HFState * hf)
{
  if (GFS_HAS_DATA (cell, h)) {
    FttCell * neighbor = ftt_cell_neighbor (cell, hf->d);
    if (!neighbor) /* boundary cell is a one-sided solid boundary: give up */
      return;
    g_assert (GFS_CELL_IS_BOUNDARY (neighbor));
    /* fixme: 
     * The boundary condition is not evaluated in the cell
     * containing the contact line.
     */
    gdouble theta = contact_angle_bc (cell, hf);
    if (theta == M_PI/2.)
      GFS_VALUE (neighbor, h) = GFS_VALUE (cell, h);
    else {
      /* fixme?: 
       * The tangential bc saturates at SLOPE_MAX. Curvatures defined
       * using parabolic interpolation are not consistent when the
       * ordinates differ too much. This is not a problem if the interface
       * is well-resolved (because the curvature will then be defined
       * using the heights in the other direction, which leads to a
       * well-defined curvature with the correct contact angle condition).
       *
       * If the interface is not well-resolved and if the desired contact
       * angle is smaller than THETA_MIN (or larger than M_PI -
       * THETA_MIN), the contact angle will saturate at THETA_MIN = atan (1./SLOPE_MAX).
       */
      gdouble cotantheta = (theta < THETA_MIN ? SLOPE_MAX : 
			    theta > M_PI - THETA_MIN ? - SLOPE_MAX :
			    1./tan(theta));
      GFS_VALUE (neighbor, h) = GFS_VALUE (cell, h) + cotantheta;
    }
  }
}

static void height_contact_tangential_bc (FttCell * cell, HFState * hf)
{
  if (is_interfacial (cell, hf->f)) {
    contact_angle_height (cell, hf->hb, hf);
    contact_angle_height (cell, hf->ht, hf);
  }
}

static void box_periodic_bc (GfsBox * box, HFState * hf)
{
  for (hf->d = 2*hf->c; hf->d <= 2*hf->c + 1; hf->d++)
    if (GFS_IS_BOUNDARY_PERIODIC (box->neighbor[hf->d]))
      ftt_cell_traverse_boundary (box->root, hf->d, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
				  (FttCellTraverseFunc) height_periodic_bc, hf);
}

static void box_contact_bc (GfsBox * box, HFState * hf)
{
  /* fixme: 2D only */
  for (hf->d = 0; hf->d < FTT_NEIGHBORS; hf->d++)
    if (GFS_IS_BOUNDARY (box->neighbor[hf->d]) && 
	!GFS_IS_BOUNDARY_PERIODIC (box->neighbor[hf->d])) {
      hf->angle = gfs_boundary_lookup_bc (GFS_BOUNDARY (box->neighbor[hf->d]), hf->f);
      if (!GFS_IS_BC_ANGLE (hf->angle))
	hf->angle = NULL; /* symmetry i.e. angle = PI/2 */
      FttCellTraverseFunc contact_bc = 
	(FttCellTraverseFunc) (hf->d/2 == hf->c ? 
			       height_contact_normal_bc : height_contact_tangential_bc);
      ftt_cell_traverse_boundary (box->root, hf->d, FTT_POST_ORDER, FTT_TRAVERSE_ALL, -1,
				  contact_bc, hf);
    }
}

static void remaining_boundary_height_undefined (FttCell * cell, HFState * hf)
{
  GfsVariable * h = boundary_hit (cell, hf);
  if (h)
    GFS_VALUE (cell, h) = GFS_NODATA;
}

static gboolean height_normal (FttCell * cell, GfsVariable * v, FttVector * m)
{
  gdouble slope = G_MAXDOUBLE;
#ifdef FTT_2D
  GfsVariableTracerVOFHeight * t = GFS_VARIABLE_TRACER_VOF_HEIGHT (v);
  FttComponent c;
  m->x = 0.;
  m->y = 1.;
  for (c = 0; c < 2; c++) {
    gdouble orientation;
    GfsVariable * hv = gfs_closest_height (cell, t, c, &orientation);
    if (hv != NULL && fabs (GFS_VALUE (cell, hv)) <= 1.) {
      gdouble H = GFS_VALUE (cell, hv);
      gdouble x[2], h[2];
      FttComponent oc = FTT_ORTHOGONAL_COMPONENT (c);
      h[0] = neighboring_column (cell, hv, c, orientation, 2*oc, &x[0]);
      if (h[0] == GFS_NODATA)
	continue;
      h[1] = neighboring_column (cell, hv, c, orientation, 2*oc + 1, &x[1]);
      if (h[1] == GFS_NODATA)
	continue;
      x[1] = - x[1];
      
      gdouble det = x[0]*x[1]*(x[0] - x[1]), a = x[1]*(h[0] - H), b = x[0]*(h[1] - H);
      gdouble hx = (x[0]*b - x[1]*a)/det;
      if (fabs (hx) < slope) {
	slope = fabs (hx);
	(&m->x)[c] = orientation;
	(&m->x)[(c + 1) % 2] = - hx;
      }
    }
  }
#else /* 3D */
  g_assert_not_implemented ();
#endif /* 3D */
  return slope < G_MAXDOUBLE;
}

static void vof_height_plane (FttCell * cell, GfsVariable * v)
{
  if (FTT_CELL_IS_LEAF (cell)) {
    GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
    gdouble f = GFS_VALUE (cell, v);
    FttComponent c;

    THRESHOLD (f);
    if (GFS_IS_FULL (f)) {
      for (c = 1; c < FTT_DIMENSION; c++)
	GFS_VALUE (cell, t->m[c]) = 0.;
      GFS_VALUE (cell, t->m[0]) = 1.;
      GFS_VALUE (cell, t->alpha) = f;
    }
    else {
      FttVector m;
      gdouble n = 0.;

      if (!height_normal (cell, v, &m))
	myc_normal (cell, v, &m);
      for (c = 0; c < FTT_DIMENSION; c++)
	n += fabs ((&m.x)[c]);
      if (n > 0.)
	for (c = 0; c < FTT_DIMENSION; c++)
	  (&m.x)[c] /= n;
      else /* fixme: this is a small fragment */
	m.x = 1.;
      for (c = 0; c < FTT_DIMENSION; c++)
	GFS_VALUE (cell, t->m[c]) = (&m.x)[c];
      GFS_VALUE (cell, t->alpha) = gfs_plane_alpha (&m, f);
    }
  }
}

static void variable_tracer_vof_height_update (GfsVariable * v, GfsDomain * domain)
{
  gfs_domain_cell_traverse (domain,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) v->fine_coarse, v);
  gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, v);
  GSList * j = GFS_VARIABLE_TRACER_VOF (v)->concentrations->items;
  while (j) {
    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, j->data);
    j = j->next;
  }

  /* update height functions */
  GfsVariableTracerVOFHeight * h = GFS_VARIABLE_TRACER_VOF_HEIGHT (v);
  HFState hf;
  hf.f = v;
  for (hf.c = 0; hf.c < FTT_DIMENSION; hf.c++) {
    hf.hb = h->hb[hf.c];
    hf.ht = h->ht[hf.c];
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) undefined_height, &hf);
    gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
					(FttCellTraverseFunc) height, &hf,
					is_interfacial, hf.f);
    
    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, hf.hb);
    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, hf.ht);
    /* update periodic boundaries first */
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_periodic_bc, &hf);
    /* apply bc again to make sure periodic bcs are in sync */
    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, hf.hb);
    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, hf.ht);
    /* apply contact angle bcs */
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_contact_bc, &hf);

    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) remaining_boundary_height_undefined, &hf);
  }

  /* update normals and alpha */
  GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
  guint l, depth = gfs_domain_depth (domain);
  FttComponent c;
  for (l = 0; l <= depth; l++) {
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, l,
			      (FttCellTraverseFunc) vof_height_plane, v);
    for (c = 0; c < FTT_DIMENSION; c++)
      gfs_domain_bc (domain, FTT_TRAVERSE_LEVEL, l, t->m[c]);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEVEL, l, t->alpha);
  }
}

static void variable_tracer_vof_height_destroy (GtsObject * o)
{
  GfsVariableTracerVOF * v = GFS_VARIABLE_TRACER_VOF (o);

  if (v->alpha) {
    FttComponent c;
    for (c = 0; c < FTT_DIMENSION; c++) {
      gts_object_destroy (GTS_OBJECT (GFS_VARIABLE_TRACER_VOF_HEIGHT (v)->hb[c]));
      gts_object_destroy (GTS_OBJECT (GFS_VARIABLE_TRACER_VOF_HEIGHT (v)->ht[c]));
    }
  }

  (* GTS_OBJECT_CLASS (gfs_variable_tracer_vof_height_class ())->parent_class->destroy) (o);
}

static void undefined_coarse_fine (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  int i;
  ftt_cell_children (parent, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i])
      GFS_VALUE (child.c[i], v) = GFS_NODATA;
}

static void variable_tracer_vof_height_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_tracer_vof_height_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsVariable * v = GFS_VARIABLE (*o);
  GfsVariableTracerVOFHeight * t = GFS_VARIABLE_TRACER_VOF_HEIGHT (v);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++) {
    static gchar index[][2] = {"x", "y", "z"};
    gchar * name = g_strdup_printf ("%s_Hb%s", v->name, index[c]);
    gchar * description = g_strdup_printf ("%s-component (bottom) of the height function for the interface defined by %s",
					   index[c], v->name);
    t->hb[c] = gfs_domain_get_or_add_variable (v->domain, name, description);
    t->hb[c]->fine_coarse = no_coarse_fine;
    t->hb[c]->coarse_fine = undefined_coarse_fine;
    g_free (name);
    g_free (description);

    name = g_strdup_printf ("%s_Ht%s", v->name, index[c]);
    description = g_strdup_printf ("%s-component (top) of the height function for the interface defined by %s",
				   index[c], v->name);
    t->ht[c] = gfs_domain_get_or_add_variable (v->domain, name, description);
    t->ht[c]->fine_coarse = no_coarse_fine;
    t->ht[c]->coarse_fine = undefined_coarse_fine;
    g_free (name);
    g_free (description);
  }
}

static void variable_tracer_vof_height_class_init (GtsObjectClass * klass)
{
  GFS_VARIABLE_TRACER_VOF_CLASS (klass)->update = variable_tracer_vof_height_update;
  klass->destroy = variable_tracer_vof_height_destroy;
  klass->read = variable_tracer_vof_height_read;
}

GfsVariableTracerVOFClass * gfs_variable_tracer_vof_height_class (void)
{
  static GfsVariableTracerVOFClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsVariableTracerVOFHeight",
      sizeof (GfsVariableTracerVOFHeight),
      sizeof (GfsVariableTracerVOFClass),
      (GtsObjectClassInitFunc) variable_tracer_vof_height_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_tracer_vof_class ()), &info);
  }

  return klass;
}

/**
 * gfs_closest_height:
 * @cell: a #FttCell.
 * @t: a #GfsVariableTracerVOFHeight.
 * @c: x, y or z.
 * @orientation: the orientation or %NULL.
 *
 * Returns: the variable containing the height value (in direction @c)
 * closest to the interface, in which case orientation is set to 1 or
 * -1, or %NULL if no heights are defined for this cell.
 */
GfsVariable * gfs_closest_height (FttCell * cell, 
				  GfsVariableTracerVOFHeight * t,
				  FttComponent c,
				  gdouble * orientation)
{
  g_return_val_if_fail (cell != NULL, NULL);
  g_return_val_if_fail (t != NULL, NULL);

  GfsVariable * hv = NULL;
  gdouble o = 0.;
  if (cell == NULL)
    return NULL;
  if (GFS_HAS_DATA (cell, t->hb[c])) {
    hv = t->hb[c]; o = 1.;
    if (GFS_HAS_DATA (cell, t->ht[c]) && 
	fabs (GFS_VALUE (cell, t->ht[c])) < fabs (GFS_VALUE (cell, t->hb[c]))) {
      hv = t->ht[c]; o = -1.;
    }
  }
  else if (GFS_HAS_DATA (cell, t->ht[c])) {
    hv = t->ht[c]; o = -1.;
  }
  if (orientation) *orientation = o;
  return hv;
}

/** \endobject{GfsVariableTracerVOFHeight} */

static gdouble interface_fractions (FttVector m, gdouble alpha, FttDirection d)
{
  gdouble f;
#if FTT_2D
  FttComponent c1 = d > 1, c2 = !c1;
  if ((&m.x)[c2] == 0) {
    gdouble sign = (d % 2 ? -1. : 1.);
    f = (sign*(&m.x)[c1] > 0.) ? 0. : 1.;
  }
  else {
    f = (alpha-(&m.x)[c1]*!(d % 2))/(&m.x)[c2];
    if(f < 0.) f = 0.; else if (f > 1.) f = 1.;
    if((&m.x)[c2] < 0.) f = 1.-f;
  }
#else /* 3D */
  FttComponent c1 = (d/2+1) % 3, c2 = (d/2+2) % 3;
  FttVector mp;
  mp.x = (&m.x)[c1];
  mp.y = (&m.x)[c2];
  f = gfs_line_area (&mp, d % 2 ? alpha : alpha - (&m.x)[d/2]);
#endif /* 3D */
  return f;
}

gdouble gfs_vof_face_fraction (const FttCellFace * face,
			       GfsVariableTracerVOF * t)
{
  g_return_val_if_fail (face != NULL, 0.);
  g_return_val_if_fail (t != NULL, 0.);

  GfsVariable * v = GFS_VARIABLE (t);
  gdouble vright, vleft = GFS_VALUE (face->cell, v);

  if (vleft == 0.)
    return 0.;
  else if (vleft != 1.0) {
    FttComponent c;
    FttVector m;
    gdouble alpha;
    for(c = 0; c < FTT_DIMENSION; c++)
      (&m.x)[c] = GFS_VALUE (face->cell, t->m[c]);
    alpha = GFS_VALUE (face->cell, t->alpha);
    vleft = interface_fractions (m, alpha, face->d);
  }

  vright = GFS_VALUE (face->neighbor, v);
  if (vright == 0.)
    return 0.;
  else if (vright != 1.0) {
    FttComponent c;
    FttVector m;
    gdouble alpha;
    for(c = 0; c < FTT_DIMENSION; c++)
      (&m.x)[c] = GFS_VALUE (face->neighbor,t->m[c]);
    alpha = GFS_VALUE (face->neighbor, t->alpha);
    if (ftt_face_type (face) == FTT_FINE_COARSE) {
      FttVector p, o, q;
      ftt_face_pos (face, &p);
      ftt_cell_pos (face->neighbor, &o);
      ftt_cell_pos (face->cell, &q);
      gdouble h = ftt_cell_size (face->neighbor);
      (&p.x)[face->d/2] += face->d % 2 ? -h/4. : h/4.;
      for (c = 0; c < FTT_DIMENSION; c++)
	alpha -= (&m.x)[c]*(0.25 + ((&p.x)[c] - (&o.x)[c])/h);
      alpha *= 2.;
    }
    vright = interface_fractions (m, alpha, FTT_OPPOSITE_DIRECTION (face->d));
  }
  return sqrt(vleft*vright);
}
