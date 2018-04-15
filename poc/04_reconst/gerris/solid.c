/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001 National Institute of Water and Atmospheric Research
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
 * \brief Solid boundaries.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "solid.h"
#include "vof.h"
#include "variable.h"

/**
 * gfs_cell_fluid:
 * @cell: a #FttCell.
 * 
 * Sets @cell and all its descendants as full fluid cells.
 */
void gfs_cell_fluid (FttCell * cell)
{
  g_return_if_fail (cell != NULL);

  if (GFS_STATE (cell)->solid) {
    g_free (GFS_STATE (cell)->solid);
    GFS_STATE (cell)->solid = NULL;
  }

  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCellChildren child;
    guint i;
 
    ftt_cell_children (cell, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i])
	gfs_cell_fluid (child.c[i]);
  }
}

typedef struct {
  GtsPoint p[4];
  GfsSegment s[4];
} CellFace;

static void face_fractions (CellFace * f, GfsSolidVector * solid, FttVector * h)
{
  static guint etod[] = { 3, 0, 2, 1 };
  guint k, m;
  gboolean ins;
  guint o = 0;
  GtsPoint r[2];
  gdouble a, x0 = f->p[0].x, y0 = f->p[0].y;
  
  solid->a = 0.;
  solid->cm.x = solid->cm.y = solid->cm.z = 0.;
  solid->ca.z = 0.;
      
  for (m = 0; m < 4 && f->s[m].n == 0; m++);
  ins = f->s[m].inside < 0;
  for (k = m; k < m + 4; k++) {
    guint i = k % 4, i1 = (i + 1) % 4;
    gdouble x1 = f->p[i].x - x0, y1 = f->p[i].y - y0, x2 = f->p[i1].x - x0, y2 = f->p[i1].y - y0;
    if (f->s[i].n > 0) {
      g_assert (ins == (f->s[i].inside < 0));
      solid->s[etod[i]] = ins ? f->s[i].x : 1. - f->s[i].x;
      r[o].x = x1 + f->s[i].x*(x2 - x1);
      r[o].y = y1 + f->s[i].x*(y2 - y1);
      if (ins) {
	x2 = r[o].x; y2 = r[o].y;
      }
      else {
	x1 = r[o].x; y1 = r[o].y;
      }
      solid->a += (x1 + x2)*(y2 - y1);
      solid->cm.x += (x1 - x2)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
      solid->cm.y += (y2 - y1)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
      o++;
      if (o == 2) {
	o = 0;
	if (ins) {
	  x1 = r[1].x; y1 = r[1].y;
	  x2 = r[0].x; y2 = r[0].y;	    
	}
	else {
	  x1 = r[0].x; y1 = r[0].y;
	  x2 = r[1].x; y2 = r[1].y;	    
	}
	solid->a += (x1 + x2)*(y2 - y1);
	solid->cm.x += (x1 - x2)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
	solid->cm.y += (y2 - y1)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
	solid->ca.x = (x1 + x2)/2.;
	solid->ca.y = (y1 + y2)/2.;
      }
      ins = !ins;
    }
    else if (ins) {
      solid->s[etod[i]] = 1.;
      solid->a += (x1 + x2)*(y2 - y1);
      solid->cm.x += (x1 - x2)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
      solid->cm.y += (y2 - y1)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
    }
    else
      solid->s[etod[i]] = 0.;
  }
  
  a = solid->a < 0. ? 0. : solid->a/(2.*h->x*h->y);
  solid->ca.x += x0;
  solid->ca.y += y0;
  if (a > 1e-4) {
    solid->cm.x = x0 + solid->cm.x/(3.*solid->a);
    solid->cm.y = y0 + solid->cm.y/(3.*solid->a);
  }
  else {
    guint n = 0;

    solid->cm.x = solid->cm.y = 0.;
    for (m = 0; m < 4 && f->s[m].n == 0; m++);
    ins = f->s[m].inside < 0;
    for (k = m; k < m + 4; k++) {
      guint i = k % 4, i1 = (i + 1) % 4;
      gdouble x1 = f->p[i].x - x0, y1 = f->p[i].y - y0, x2 = f->p[i1].x - x0, y2 = f->p[i1].y - y0;
      if (f->s[i].n > 0) {
	gdouble x = x1 + f->s[i].x*(x2 - x1);
	gdouble y = y1 + f->s[i].x*(y2 - y1);

	g_assert (ins == (f->s[i].inside < 0));
	solid->cm.x += x;
	solid->cm.y += y;
	n++;
	if (ins) {
	  solid->cm.x += x1;
	  solid->cm.y += y1;
	  n++;
	}
	ins = !ins;
      }
      else if (ins) {
	solid->cm.x += x1;
	solid->cm.y += y1;
	n++;
      }
    }
    g_assert (n > 0);
    solid->cm.x = x0 + solid->cm.x/n;
    solid->cm.y = y0 + solid->cm.y/n;
  }
  solid->a = a;
}

static void face_new (CellFace * f, FttCell * cell, GfsGenericSurface * s, FttVector * h)
{
  FttVector p;
  guint i;

  ftt_cell_pos (cell, &p);
  f->p[0].x = p.x - h->x/2.; f->p[0].y = p.y - h->y/2.; f->p[0].z = 0.;
  f->p[1].x = p.x + h->x/2.; f->p[1].y = p.y - h->y/2.; f->p[1].z = 0.;
  f->p[2].x = p.x + h->x/2.; f->p[2].y = p.y + h->y/2.; f->p[2].z = 0.;
  f->p[3].x = p.x - h->x/2.; f->p[3].y = p.y + h->y/2.; f->p[3].z = 0.;

  for (i = 0; i < 4; i++) {
    f->s[i].E = &f->p[i];
    f->s[i].D = &f->p[(i + 1) % 4];
    gfs_surface_segment_intersection (s, cell, &f->s[i]);
  }
}

static gboolean solid_face_is_thin (CellFace * f)
{
  guint odd = 0, even = 0, i;

  for (i = 0; i < 4; i++)
    if (f->s[i].n) {
      if (f->s[i].n % 2 != 0)
	odd++;
      else
	even++;
    }
  if (odd == 2 && even == 1) {
    for (i = 0; i < 4; i++)
      if (f->s[i].n % 2 != 0 && f->s[(i + 2) % 4].n % 2 != 0)
	return FALSE;
    return TRUE;
  }
  return (odd > 2 || even > 1);
}

/**
 * gfs_set_2D_solid_fractions_from_surface:
 * @cell: a #FttCell.
 * @s: a #GfsGenericSurface.
 *
 * Sets the 2D volume fractions of @cell cut by @s.
 *
 * Returns: %TRUE if the cell is thin, %FALSE otherwise;
 */
gboolean gfs_set_2D_solid_fractions_from_surface (FttCell * cell,
						  GfsGenericSurface * s)
{
  GfsSolidVector * solid;
  FttVector h;
  CellFace f;
  guint i, n1 = 0;
  gboolean thin = FALSE;

  g_return_val_if_fail (cell != NULL, FALSE);
  g_return_val_if_fail (s != NULL, FALSE);

  h.x = h.y = ftt_cell_size (cell);
  face_new (&f, cell, s, &h);
  
  for (i = 0; i < 4; i++)
    if (f.s[i].n % 2 != 0) {
      f.s[i].x /= f.s[i].n;
      n1++;
    }
    else
      f.s[i].n = 0;

  solid = GFS_STATE (cell)->solid;
  switch (n1) {
  case 0:
    break;
  case 4:
    thin = TRUE;
    /* fall through */
  case 2: {
    if (!solid)
      GFS_STATE (cell)->solid = solid = g_malloc0 (sizeof (GfsSolidVector));
    face_fractions (&f, solid, &h);
    if (solid->a == 1.) {
      g_free (solid);
      GFS_STATE (cell)->solid = NULL;
    }
    break;
  }
  default: {
    FttVector p;
    ftt_cell_pos (cell, &p);
    g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	   "the surface may not be closed (n1 = %d)\n"
	   "at (%g,%g,%g)", n1, p.x, p.y, p.z);
  }
  }
  return thin;
}

typedef struct {
  gboolean destroy_solid;
  FttCellCleanupFunc cleanup;
  gpointer data;
  GfsVariable * status;
  guint thin;
  GSList * solid_boxes;
} InitSolidParams;

static gboolean thin_cell_is_solid (FttCell * cell)
{
  gdouble sum = 0.;
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++)
    sum += GFS_STATE (cell)->solid->s[d];
  return (sum < FTT_NEIGHBORS/2);
}

static void deal_with_thin_cell (FttCell * cell, InitSolidParams * p)
{
  cell->flags |= GFS_FLAG_THIN;
  if (thin_cell_is_solid (cell))
    GFS_VALUE (cell, p->status) = GFS_STATUS_SOLID;
  else {
    GfsSolidVector * solid = GFS_STATE (cell)->solid;
    FttDirection d;
    for (d = 0; d < FTT_NEIGHBORS; d++)
      solid->s[d] = (solid->s[d] > 0.5);
    solid->a = 1.;
    ftt_cell_pos (cell, &solid->cm);
    solid->ca = solid->cm;
  }
}

#if FTT_2D /* 2D */

static void set_solid_fractions_from_surface (FttCell * cell,
					      GfsGenericSurface * s,
					      InitSolidParams * p)
{
  if (gfs_set_2D_solid_fractions_from_surface (cell, s)) {
    p->thin++;
    deal_with_thin_cell (cell, p);
  }
  else if (GFS_STATE (cell)->solid && GFS_STATE (cell)->solid->a == 0.)
    GFS_VALUE (cell, p->status) = GFS_STATUS_SOLID;
}

/**
 * gfs_solid_is_thin:
 * @cell: a #FttCell.
 * @s: a #GfsGenericSurface.
 *
 * @s is "thin" relative to @cell if the miminum distance between
 * non-connected faces of @s cutting @cell is smaller than the size of
 * @cell (see doc/figures/thin.fig).
 *
 * Returns: %TRUE if @s is a thin surface, %FALSE otherwise.
 */
gboolean gfs_solid_is_thin (FttCell * cell, GfsGenericSurface * s)
{
  CellFace f;
  FttVector h;

  g_return_val_if_fail (cell != NULL, FALSE);
  g_return_val_if_fail (s != NULL, FALSE);

  h.x = h.y = ftt_cell_size (cell);
  face_new (&f, cell, s, &h);
  return solid_face_is_thin (&f);
}

#else /* 3D */
#include "isocube.h"

typedef struct {
  GtsPoint p[8];
  GfsSegment s[12];
} CellCube;

static void rotate (CellFace * f, FttVector * h, FttComponent c)
{
  guint i;

  switch (c) {
  case FTT_X: 
    for (i = 0; i < 4; i++) {
      f->p[i].x = f->p[i].y; f->p[i].y = f->p[i].z;
    }
    h->x = h->y; h->y = h->z;
    break;
  case FTT_Y:
    for (i = 0; i < 4; i++)
      f->p[i].y = f->p[i].z;
    h->y = h->z;
    break;
  case FTT_Z:
    break;
  default:
    g_assert_not_reached ();
  }
}

static void cell_size (FttCell * cell, FttVector * h)
{
  h->x = h->y = ftt_cell_size (cell);
  h->z = h->x;
}

/* Returns: the number of closed loops for the given isocube
 * 
 * Fixme: this algorithm is not correct in general. This has no
 * consequence for this particular application because we also check
 * that the isosurface is "planar" together with use of topology() in
 * set_solid_fractions_from_surface3D(), however this is important in
 * general.
 *
 * The bug is triggered for certain configurations of "non-planar"
 * isosurfaces.
 */
static guint topology (CellCube * cube)
{
  guint l, nl = 0;
  gboolean used[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
  
  for (l = 0; l < 12; l++) {
    guint nv = 0, e = l, cut = cube->s[e].n % 2;
    
    while (cut && !used[e]) {
      guint m = 0, * ne = connect[e][cube->s[e].inside > 0];

      nv++;
      used[e] = TRUE;
      cut = 0;
      while (m < 3 && !cut) {
	e = ne[m++];
	cut = cube->s[e].n % 2;
      }
    }
    if (nv > 2)
      nl++;
  }
  return nl;
}

static void cube_new (CellCube * cube, FttCell * cell, GfsGenericSurface * s, FttVector * o, FttVector * h)
{
  guint i;

  for (i = 0; i < FTT_DIMENSION; i++)
    (&o->x)[i] -= (&h->x)[i]/2.;
  for (i = 0; i < 8; i++) { /* for each vertex of the cube */
    cube->p[i].x = o->x + h->x*vertex[i].x;
    cube->p[i].y = o->y + h->y*vertex[i].y;
    cube->p[i].z = o->z + h->z*vertex[i].z;
  }

  for (i = 0; i < 12; i++) {
    cube->s[i].E = &cube->p[edge1[i][0]];
    cube->s[i].D = &cube->p[edge1[i][1]];
    gfs_surface_segment_intersection (s, cell, &cube->s[i]);
  }
}

static void set_solid_fractions_from_surface (FttCell * cell, 
					      GfsGenericSurface * surface, 
					      InitSolidParams * p)
{
  GfsSolidVector * solid = GFS_STATE (cell)->solid;
  CellCube cube;
  FttVector o, ca = {0., 0., 0.}, h;
  guint i, n1 = 0;
  gint inside[8] = {0,0,0,0,0,0,0,0};
  gboolean planar = TRUE;

  ftt_cell_pos (cell, &o);
  cell_size (cell, &h);
  cube_new (&cube, cell, surface, &o, &h);

  for (i = 0; i < 12; i++) { /* for each edge of the cube */
    GfsSegment * s = &cube.s[i];
    if (cube.s[i].n % 2 != 0) { /* only for odd number of intersections */
      guint j = edge1[i][0], k = edge1[i][1];

      /* intersection vertex position is the average of all the n[i] intersections */
      s->x /= s->n;

      /* average of all intersections */
      ca.x += (1. - s->x)*cube.p[j].x + s->x*cube.p[k].x;
      ca.y += (1. - s->x)*cube.p[j].y + s->x*cube.p[k].y;
      ca.z += (1. - s->x)*cube.p[j].z + s->x*cube.p[k].z;

      g_assert (inside[j] == 0 || inside[j] == s->inside);
      g_assert (inside[k] == 0 || inside[k] == - s->inside);
      inside[j] = s->inside;
      inside[k] = - s->inside;
      n1++;
    }
    else
      s->n = 0;
  }

  if (n1 == 0) { /* no intersections */
    if (solid) {
      g_free (solid);
      GFS_STATE (cell)->solid = NULL;      
    }
    return;
  }

  if (!solid)
    GFS_STATE (cell)->solid = solid = g_malloc0 (sizeof (GfsSolidVector));

  /* compute face fractions */
  for (i = 0; i < FTT_NEIGHBORS; i++) {
    CellFace f;
    guint j, n2;

    n2 = 0;
    for (j = 0; j < 4; j++) { /* initialise face i */
      GfsSegment * s = &cube.s[face[i][j][0]];

      f.p[j] = cube.p[face_v[i][j]];
      f.s[j].n = s->n;
      if (f.s[j].n) n2++;
      if (face[i][j][1]) {
	f.s[j].x = 1. - s->x;
	f.s[j].inside = - s->inside;
      }
      else {
	f.s[j].x = s->x;
	f.s[j].inside = s->inside;
      }
    }

    switch (n2) {
    case 0: { /* the face is not cut */
      gint ins = 0;

      /* checks whether the face vertices are inside or outside */
      for (j = 0; j < 4; j++) {
	gint k = inside[face_v[i][j]];
	if (k) {
	  g_assert (ins == 0 || ins == k);
	  ins = k;
	}
      }
      g_assert (ins != 0);
      solid->s[i] = ins > 0 ? 0. : 1.;
      break;
    }
    case 4:
      planar = FALSE;
      /* fall through */
    case 2: { /* the face is cut 2 or 4 times */
      GfsSolidVector sol;
      FttVector h1;

      h1 = h;
      rotate (&f, &h1, i/2);
      face_fractions (&f, &sol, &h1);
      solid->s[i] = sol.a;
      break;
    }
    default: {
      FttVector p;
      ftt_cell_pos (cell, &p);
      g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	     "the surface may not be closed (n2 = %d)\n"
	     "at (%g,%g,%g)", n2, p.x, p.y, p.z);
    }
    }
  }

  /* now compute cell fraction, center of area, center of mass */
  ca.x /= n1; ca.y /= n1; ca.z /= n1; 
  solid->ca = ca;
  if (planar && topology (&cube) == 1) {
    FttVector m;
    gdouble alpha, n = 0.;
    gboolean sym[FTT_DIMENSION];
    FttComponent c;

    for (c = 0; c < FTT_DIMENSION; c++) {
      (&ca.x)[c] = ((&ca.x)[c] - (&o.x)[c])/(&h.x)[c];
      (&m.x)[c] = solid->s[2*c + 1] - solid->s[2*c];
      if ((&m.x)[c] < 0.) {
	(&m.x)[c] = - (&m.x)[c];
	(&ca.x)[c] = 1. - (&ca.x)[c];
	sym[c] = TRUE;
      }
      else
	sym[c] = FALSE;
      n += (&m.x)[c];
    }
    if (n == 0.) { /* this is a fluid or solid cell */
      for (c = 1; c < FTT_NEIGHBORS; c++)
	g_assert (solid->s[c] == solid->s[0]);
      if (solid->s[0] == 1.) { /* fluid */
	g_free (solid);
	GFS_STATE (cell)->solid = NULL;
	return;
      }
      else { /* solid */
	solid->a = 0.;
	solid->cm.x = solid->cm.y = solid->cm.z = 0.;
      }
    }
    else {
      m.x /= n; m.y /= n; m.z /= n;
      alpha = m.x*ca.x + m.y*ca.y + m.z*ca.z;
      solid->a = gfs_plane_volume (&m, alpha);
      gfs_plane_center (&m, alpha, solid->a, &solid->cm);
    }
    for (c = 0; c < FTT_DIMENSION; c++)
      (&solid->cm.x)[c] = (&o.x)[c] + 
	(sym[c] ? 1. - (&solid->cm.x)[c] : (&solid->cm.x)[c])*(&h.x)[c];
  }
  else { /* this is a "thin" cell */
    p->thin++;
    deal_with_thin_cell (cell, p);
  }
  if (solid->a == 0.)
    GFS_VALUE (cell, p->status) = GFS_STATUS_SOLID;
}

/**
 * gfs_solid_is_thin:
 * @cell: a #FttCell.
 * @s: a #GfsGenericSurface.
 *
 * @s is "thin" relative to @cell if the miminum distance between
 * non-connected faces of @s cutting @cell is smaller than the size of
 * @cell (see doc/figures/thin.fig).
 *
 * Returns: %TRUE if @s is a thin surface, %FALSE otherwise.
 */
gboolean gfs_solid_is_thin (FttCell * cell, GfsGenericSurface * s)
{
  CellCube cube;
  FttVector o, h;
  guint i;

  g_return_val_if_fail (cell != NULL, FALSE);
  g_return_val_if_fail (s != NULL, FALSE);

  ftt_cell_pos (cell, &o);
  cell_size (cell, &h);
  cube_new (&cube, cell, s, &o, &h);
  for (i = 0; i < FTT_NEIGHBORS; i++) {
    CellFace f;
    guint j;

    for (j = 0; j < 4; j++)
      f.s[j].n = cube.s[face[i][j][0]].n;
    if (solid_face_is_thin (&f))
      return TRUE;
  }
  return (topology (&cube) > 1);
}

#endif /* 3D */

static gdouble solid_sa (GfsSolidVector * s)
{
  gdouble sa2 = 0.;
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++) {
    gdouble n = s->s[2*c] - s->s[2*c + 1];

    sa2 += n*n;
  }
  return sqrt (sa2);
}

/**
 * gfs_cell_init_solid_fractions_from_children:
 * @cell: a #FttCell.
 *
 * Uses the values of the solid fractions of the children of @cell to
 * compute the values of its solid fractions.
 *
 * This function fails if @cell is a leaf of the cell tree.  
 */
void gfs_cell_init_solid_fractions_from_children (FttCell * cell)
{
  FttCellChildren child;
  guint i, j;
  gdouble w = 0., wa = 0.;
  gboolean cell_is_solid = TRUE;
  gboolean cell_is_mixed = FALSE;
  FttVector cm = { 0., 0., 0.};
  FttVector ca = { 0., 0., 0.};

  g_return_if_fail (cell != NULL);
  g_return_if_fail (!FTT_CELL_IS_LEAF (cell));

  ftt_cell_children (cell, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i]) {
      if (GFS_IS_FLUID (child.c[i])) {
	FttVector p;

	w += 1.;
	ftt_cell_pos (child.c[i], &p);
	cm.x += p.x; cm.y += p.y; cm.z += p.z;
	cell_is_solid = FALSE;
      }
      else {
	GfsSolidVector * solid = GFS_STATE (child.c[i])->solid;
	gdouble sa = solid_sa (solid) + 1e-9;

	w += solid->a; wa += sa;
	cm.x += solid->cm.x*solid->a;
	cm.y += solid->cm.y*solid->a;
	cm.z += solid->cm.z*solid->a;
	ca.x += solid->ca.x*sa;
	ca.y += solid->ca.y*sa;
	ca.z += solid->ca.z*sa;
	cell_is_mixed = TRUE;
      }
    }

  if (cell_is_mixed) {
    GfsSolidVector * solid = GFS_STATE (cell)->solid;

    if (solid == NULL)
      GFS_STATE (cell)->solid = solid = g_malloc0 (sizeof (GfsSolidVector));

    solid->a = w/FTT_CELLS;
    g_assert (wa > 0.);
    solid->ca.x = ca.x/wa;
    solid->ca.y = ca.y/wa;
    solid->ca.z = ca.z/wa;
    if (w > 0.) {
      solid->cm.x = cm.x/w;
      solid->cm.y = cm.y/w;
      solid->cm.z = cm.z/w;
    }
    else
      ftt_cell_pos (cell, &solid->cm);

    for (i = 0; i < FTT_NEIGHBORS; i++) {
      guint n = ftt_cell_children_direction (cell, i, &child);

      w = 0.;
      for (j = 0; j < n; j++)
	if (child.c[j])
	  w += GFS_IS_FLUID (child.c[j]) ? 1. : GFS_STATE (child.c[j])->solid->s[i];
      solid->s[i] = w/n;
    }
  }
  else { /* !cell_is_mixed */
    if (GFS_STATE (cell)->solid) {
      g_free (GFS_STATE (cell)->solid);
      GFS_STATE (cell)->solid = NULL;
    }
    g_assert (!cell_is_solid);
  }
}

static void push_leaf (GtsFifo * fifo, FttCell * cell, FttDirection d, gdouble a,
		       GfsVariable * status)
{
  if (FTT_CELL_IS_LEAF (cell)) {
    if (!GFS_IS_MIXED (cell) && GFS_VALUE (cell, status) == GFS_STATUS_UNDEFINED) {
      GFS_VALUE (cell, status) = a;
      gts_fifo_push (fifo, cell);
    }
  }
  else {
    FttCellChildren child;
    guint i, n;
    
    n = ftt_cell_children_direction (cell, FTT_OPPOSITE_DIRECTION (d), &child);
    for (i = 0; i < n; i++)
      if (child.c[i] && !GFS_IS_MIXED (child.c[i]) && 
	  GFS_VALUE (child.c[i], status) == GFS_STATUS_UNDEFINED) {
	g_assert (FTT_CELL_IS_LEAF (child.c[i]));
	GFS_VALUE (child.c[i], status) = a;
	gts_fifo_push (fifo, child.c[i]);
      }
  }
}

static void paint_leaf (GtsFifo * fifo, gdouble a, GfsVariable * status)
{
  FttCell * cell;

  while ((cell = gts_fifo_pop (fifo))) {
    FttDirection i;
    FttCellNeighbors n;
    
    ftt_cell_neighbors (cell, &n);
    for (i = 0; i < FTT_NEIGHBORS; i++)
      if (n.c[i] && !GFS_CELL_IS_BOUNDARY (n.c[i]))
	push_leaf (fifo, n.c[i], i, a, status);
  }
}

static void paint_mixed_leaf (FttCell * cell, GfsVariable * status)
{
  if (GFS_IS_MIXED (cell)) {
    GfsSolidVector * solid = GFS_STATE (cell)->solid;
    GtsFifo * fifo;
    FttCell * n;
    FttDirection i;

    fifo = gts_fifo_new ();
    for (i = 0; i < FTT_NEIGHBORS; i++)
      if ((n = ftt_cell_neighbor (cell, i)) && !GFS_CELL_IS_BOUNDARY (n)) {
	if (solid->s[i] == 0. || solid->s[i] == 1.) {
	  push_leaf (fifo, n, i, solid->s[i] + 1., status);
	  paint_leaf (fifo, solid->s[i] + 1., status);
	}
	else if (!FTT_CELL_IS_LEAF (n)) {
	  FttCellChildren child;
	  guint j, k;
	  gdouble w = 0.;

	  k = ftt_cell_children_direction (n, FTT_OPPOSITE_DIRECTION (i), &child);
	  for (j = 0; j < k; j++)
	    if (child.c[j])
	      w += GFS_IS_FLUID (child.c[j]) ? 1. : 
		GFS_STATE (child.c[j])->solid->s[FTT_OPPOSITE_DIRECTION (i)];
	  if (w/k <= 0. || w/k >= 1.)
	    g_warning ("file %s: line %d (%s): w/k=%g solid->s[%d]=%g",
		       __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
		       w/k, i, solid->s[i]);
	  solid->s[i] = w/k;
	}
      }
    gts_fifo_destroy (fifo);
  }
}

static void solid_fractions_from_children (FttCell * cell, InitSolidParams * p)
{
  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCellChildren child;
    guint i;
    
    ftt_cell_children (cell, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i])
	solid_fractions_from_children (child.c[i], p);
    if (FTT_CELL_IS_LEAF (cell))
      /* all the children have been destroyed i.e. the cell is solid */
      GFS_VALUE (cell, p->status) = GFS_STATUS_SOLID;
    else {
      gfs_cell_init_solid_fractions_from_children (cell);
      GFS_VALUE (cell, p->status) = GFS_STATUS_UNDEFINED;
      if (!p->destroy_solid && !GFS_IS_MIXED (cell)) {
	ftt_cell_children (cell, &child);
	for (i = 0; i < FTT_CELLS; i++)
	  if (child.c[i]) {
	    if (GFS_VALUE (cell, p->status) == GFS_STATUS_UNDEFINED)
	      GFS_VALUE (cell, p->status) = GFS_VALUE (child.c[i], p->status);
	    else
	      g_assert (GFS_VALUE (cell, p->status) == GFS_VALUE (child.c[i], p->status));
	  }
      }
    }
  }
  if (p->destroy_solid && 
      GFS_VALUE (cell, p->status) == GFS_STATUS_SOLID && 
      !FTT_CELL_IS_ROOT (cell))
    ftt_cell_destroy (cell, p->cleanup, p->data);
}

static void foreach_box (GfsBox * box, InitSolidParams * p)
{
  solid_fractions_from_children (box->root, p);
  if (p->destroy_solid && GFS_VALUE (box->root, p->status) == GFS_STATUS_SOLID)
    p->solid_boxes = g_slist_prepend (p->solid_boxes, box);
}

static void match_fractions (FttCell * cell, GfsVariable * status)
{
  if (GFS_IS_MIXED (cell)) {
    FttCellNeighbors neighbor;
    GfsSolidVector * solid = GFS_STATE (cell)->solid;
    FttDirection d;

    ftt_cell_neighbors (cell, &neighbor);
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (neighbor.c[d] && !GFS_CELL_IS_BOUNDARY (neighbor.c[d])) {
	if (!FTT_CELL_IS_LEAF (neighbor.c[d])) {
	  FttCellChildren child;
	  FttDirection od = FTT_OPPOSITE_DIRECTION (d);
	  guint i, n = ftt_cell_children_direction (neighbor.c[d], od, &child);
	  gdouble s = 0.;

	  g_assert (GFS_VALUE (neighbor.c[d], status) != 1.);
	  for (i = 0; i < n; i++)
	    if (child.c[i] && GFS_VALUE (child.c[i], status) != 1.)
	      s += GFS_IS_MIXED (child.c[i]) ? GFS_STATE (child.c[i])->solid->s[od] : 1.;
	  solid->s[d] = s/n;
	}
	else if (GFS_VALUE (neighbor.c[d], status) != GFS_STATUS_SOLID) {
	  if (!GFS_IS_MIXED (neighbor.c[d]))
	    solid->s[d] = 1.;
	  else if (neighbor.c[d]->flags & GFS_FLAG_THIN)
	    solid->s[d] = GFS_STATE (neighbor.c[d])->solid->s[FTT_OPPOSITE_DIRECTION (d)];
	}
	else /* neighbor.c[d] is a solid cell */
	  solid->s[d] = 0.;
      }
  }
}

/**
 * gfs_init_solid_fractions_leaves:
 * @domain: a #GfsDomain.
 * @i: a list of #GfsSolids.
 * @status: a temporary variable or %NULL.
 *
 * Initializes the solid fractions of the leaf cells of @domain.
 *
 * Returns: the number of thin cells.
 */
guint gfs_init_solid_fractions_leaves (GfsDomain * domain,
				       GSList * i,
				       GfsVariable * status)
{
  InitSolidParams p;

  g_return_val_if_fail (domain != NULL, 0);

  p.status = status ? status : gfs_temporary_variable (domain);
  p.thin = 0;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, p.status);
  while (i) {
    gfs_domain_traverse_cut (domain, GFS_SOLID (i->data)->s, 
			     FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			     (FttCellTraverseCutFunc) set_solid_fractions_from_surface, &p);
    i = i->next;
  }
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) paint_mixed_leaf, p.status);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) match_fractions, p.status);
  if (status == NULL)
    gts_object_destroy (GTS_OBJECT (p.status));

  return p.thin;
}

/**
 * gfs_init_solid_fractions_from_children:
 * @domain: a #GfsDomain.
 * @destroy_solid: controls what to do with solid cells.
 * @cleanup: a #FttCellCleanupFunc or %NULL.
 * @data: user data to pass to @cleanup.
 * @status: the status variable.
 *
 * Initializes the solid fractions of the non-leaf cells of @domain
 * using the values of the leaf cells.
 *
 * If @destroy_solid is set to %TRUE, the cells entirely contained in
 * the solid are destroyed using @cleanup as cleanup function.  
 */
void gfs_init_solid_fractions_from_children (GfsDomain * domain,
					     gboolean destroy_solid,
					     FttCellCleanupFunc cleanup,
					     gpointer data,
					     GfsVariable * status)
{
  InitSolidParams p;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (status != NULL);

  p.destroy_solid = destroy_solid;
  p.cleanup = cleanup;
  p.data = data;
  p.status = status;
  p.solid_boxes = NULL;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) foreach_box, &p);
  g_slist_foreach (p.solid_boxes, (GFunc) gts_object_destroy, NULL);
  g_slist_free (p.solid_boxes);
  if (p.solid_boxes) {
    gfs_locate_array_destroy (domain->array);
    domain->array = gfs_locate_array_new (domain);
  }
}

/**
 * gfs_domain_init_solid_fractions:
 * @domain: a #GfsDomain.
 * @i: a list of #GfsSolids.
 * @destroy_solid: controls what to do with solid cells.
 * @cleanup: a #FttCellCleanupFunc or %NULL.
 * @data: user data to pass to @cleanup.
 * @status: a temporary variable or %NULL.
 *
 * Initializes the solid fractions of all the cells of @domain.
 *
 * If @destroy_solid is set to %TRUE, the cells entirely contained in
 * the solid are destroyed using @cleanup as cleanup function.  
 *
 * Returns: the number of thin cells.
 */
guint gfs_domain_init_solid_fractions (GfsDomain * domain,
				       GSList * i,
				       gboolean destroy_solid,
				       FttCellCleanupFunc cleanup,
				       gpointer data,
				       GfsVariable * status)
{
  GfsVariable * status1;

  g_return_val_if_fail (domain != NULL, 0);

  status1 = status ? status : gfs_temporary_variable (domain);
  guint thin = gfs_init_solid_fractions_leaves (domain, i, status1);
  gfs_init_solid_fractions_from_children (domain, destroy_solid, cleanup, data, status1);
  if (status == NULL)
    gts_object_destroy (GTS_OBJECT (status1));

  return thin;
}

static gboolean check_area_fractions (const FttCell * root)
{
  guint i, level;
  FttCellNeighbors neighbor;
  gboolean ret = TRUE;
  GfsSolidVector * solid;

  level = ftt_cell_level (root);
  ftt_cell_neighbors (root, &neighbor);
  solid = GFS_STATE (root)->solid;

  if (solid) {
    GtsBBox bb;

    ftt_cell_bbox (root, &bb);
    if (!gts_bbox_point_is_inside (&bb, &solid->cm)) {
      g_warning ("file %s: line %d (%s): cm (%g,%g,%g)/%d is not inside cell [(%g,%g,%g),(%g,%g,%g)]",
		 __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
		 solid->cm.x, solid->cm.y, solid->cm.z, ftt_cell_level (root),
		 bb.x1, bb.y1, bb.z1, 
		 bb.x2, bb.y2, bb.z2);
      ret = FALSE;
      g_assert_not_reached ();
    }
    if (!gts_bbox_point_is_inside (&bb, &solid->ca)) {
      g_warning ("file %s: line %d (%s): ca (%g,%g,%g)/%d is not inside cell [(%g,%g,%g),(%g,%g,%g)]",
		 __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
		 solid->ca.x, solid->ca.y, solid->ca.z, ftt_cell_level (root),
		 bb.x1, bb.y1, bb.z1, 
		 bb.x2, bb.y2, bb.z2);
      ret = FALSE;
      g_assert_not_reached ();
    }
  }

  for (i = 0; i < FTT_NEIGHBORS; i++)
    if (neighbor.c[i]) {
      GfsSolidVector * nsolid = GFS_STATE (neighbor.c[i])->solid;
      FttDirection oi = FTT_OPPOSITE_DIRECTION (i);

      if (ftt_cell_level (neighbor.c[i]) == level) {
	if (GFS_IS_FLUID (root)) {
	  if (!GFS_IS_FLUID (neighbor.c[i])) {
	    if (1. - nsolid->s[oi] >= 1e-10) {
	      FttVector p;
	      ftt_cell_pos (root, &p);
	      g_warning ("file %s: line %d (%s): (%g,%g,%g)/%d: s[%d]: %g",
			 __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
			 p.x, p.y, p.z, ftt_cell_level (root),
			 oi, nsolid->s[oi]);
	      ret = FALSE;
	    }
	    nsolid->s[oi] = 1.;
	  }
	}
	else if (GFS_IS_MIXED (neighbor.c[i])) {
	  if (fabs (solid->s[i] - nsolid->s[oi]) >= 1e-10) {
	    FttVector p;
	    ftt_cell_pos (root, &p);
	    g_warning ("file %s: line %d (%s): (%g,%g,%g)/%d: s[%d]: %g neighbor->s[%d]: %g",
		       __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
		       p.x, p.y, p.z, ftt_cell_level (root),
		       i, solid->s[i],
		       oi, nsolid->s[oi]);
	    ret = FALSE;
	  }
	  nsolid->s[oi] = solid->s[i];
	}
	else {
	  if (1. - solid->s[i] >= 1e-10) {
	    FttVector p;
	    ftt_cell_pos (root, &p);
	    g_warning ("file %s: line %d (%s): (%g,%g,%g)/%d: s[%d]: %g",
		       __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
		       p.x, p.y, p.z, ftt_cell_level (root),
		       i, solid->s[i]);
	    ret = FALSE;
	  }
	  solid->s[i] = 1.;
	}
      }
      else { /* fine/coarse boundary */
	g_assert (ftt_cell_level (neighbor.c[i]) == level - 1);
	if (GFS_IS_FLUID (neighbor.c[i])) {
	  if (GFS_IS_MIXED (root)) {
	    if (1. - solid->s[i] >= 1e-10) {
	      FttVector p;
	      ftt_cell_pos (root, &p);
	      g_warning ("file %s: line %d (%s): (%g,%g,%g)/%d: s[%d]: %g",
			 __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
			 p.x, p.y, p.z, ftt_cell_level (root),
			 i, solid->s[i]);
	      ret = FALSE;
	    }
	    solid->s[i] = 1.;
	  }
	}
	else if (nsolid->s[oi] == 0.) {
	  g_assert (GFS_IS_MIXED (root));
	  if (solid->s[i] >= 1e-10) {
	    FttVector p;
	    ftt_cell_pos (root, &p);
	    g_warning ("file %s: line %d (%s): (%g,%g,%g)/%d: s[%d]: %g",
		       __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
		       p.x, p.y, p.z, ftt_cell_level (root),
		       i, solid->s[i]);
	    ret = FALSE;
	  }
	  solid->s[i] = 0.;
	}
      }
    }
  
  if (!FTT_CELL_IS_LEAF (root)) {
    FttCellChildren child;

    ftt_cell_children (root, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i] && !check_area_fractions (child.c[i]))
	ret = FALSE;
  }

  return ret;
}

static void check_solid_fractions (FttCell * cell, gboolean * ret)
{
  FttCellChildren children;
  guint n;

  ftt_cell_children (cell, &children);
  if (!GFS_IS_MIXED (cell)) {
    for (n = 0; n < FTT_CELLS; n++)
      if (children.c[n] && GFS_IS_MIXED (children.c[n])) {
	g_warning ("file %s: line %d (%s): children[%d] is mixed (%g)"
		   " parent is not",
                   __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
		   n, GFS_STATE (children.c[n])->solid->a);
	*ret = FALSE;
      }
  }
  else {
    gdouble a = 0.;

    for (n = 0; n < FTT_CELLS; n++)
      if (children.c[n]) {
	if (GFS_IS_MIXED (children.c[n]))
	  a += GFS_STATE (children.c[n])->solid->a;
	else
	  a += 1.;
      }
    a /= FTT_CELLS;
    if (fabs (GFS_STATE (cell)->solid->a - a) >= 1e-10) {
      g_warning ("file %s: line %d (%s): children->a: %g parent->a: %g",
		 __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
		 a, GFS_STATE (cell)->solid->a);
	*ret = FALSE;
    }
  }
}

/**
 * gfs_cell_check_solid_fractions:
 * @root: the root #FttCell of the cell tree to check.
 * 
 * Checks the consistency of the solid fractions of each cell of the
 * cell tree relative to the neighboring solid fractions.
 *
 * Returns: %TRUE if the solid fractions are consistent, %FALSE otherwise.
 */
gboolean gfs_cell_check_solid_fractions (FttCell * root)
{
  gboolean ret = TRUE;

  g_return_val_if_fail (root != NULL, FALSE);

  ftt_cell_traverse (root, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
		     (FttCellTraverseFunc) check_solid_fractions, &ret);
  return ret & check_area_fractions (root);
}

static void save_solid (FttCell * cell, GfsVariable * c)
{
  GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, c)) = GFS_STATE (cell)->solid;
  GFS_STATE (cell)->solid = NULL;
}

static void restore_solid (FttCell * cell, gpointer * data)
{
  GfsVariable * status = data[0];
  GfsVariable * c = data[1];
  GfsSolidVector * solid = GFS_STATE (cell)->solid;

  GFS_STATE (cell)->solid = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, c));
  if (solid) {
    GFS_VALUE (cell, c) = solid->a;
    g_free (solid);
  }
  else {
    g_assert (GFS_VALUE (cell, status) == GFS_STATUS_SOLID || 
	      GFS_VALUE (cell, status) == GFS_STATUS_FLUID);
    GFS_VALUE (cell, c) = GFS_VALUE (cell, status) - 1.;
  }
}

static void set_status (FttCell * cell, gpointer * data)
{
  GfsVariable * status = data[0];
  gdouble * val = data[2];
  GFS_VALUE (cell, status) = *val;
}

static void check_status (GfsBox * box, gpointer * data)
{
  GfsVariable * status = data[0];
  if (!GFS_IS_MIXED (box->root) && GFS_VALUE (box->root, status) == GFS_STATUS_UNDEFINED) {
    GfsGenericSurface * s = data[1];
    FttVector pos;
    ftt_cell_pos (box->root, &pos);
    gdouble val = gfs_surface_point_is_inside (s, &pos) > 0 ?
      GFS_STATUS_FLUID : GFS_STATUS_SOLID;
    data[2] = &val;
    ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
		       (FttCellTraverseFunc) set_status, data);
  }
}

/**
 * gfs_domain_init_fraction:
 * @domain: a #GfsDomain.
 * @s: a surface defining the interface boundary.
 * @c: a #GfsVariable.
 *
 * Initializes the fraction @c of the interface @s contained in all
 * the cells of @domain.
 */
void gfs_domain_init_fraction (GfsDomain * domain,
			       GfsGenericSurface * s,
			       GfsVariable * c)
{
  gpointer data[3];
  GfsVariable * status;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (s != NULL);
  g_return_if_fail (c != NULL);

  status = gfs_temporary_variable (domain);

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) save_solid, c);
  GfsSolid tmp;
  tmp.s = s;
  GSList * l = g_slist_prepend (NULL, &tmp);
  gfs_domain_init_solid_fractions (domain, l, FALSE, NULL, NULL, status);
  g_slist_free (l);
  data[0] = status;
  data[1] = s;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) check_status, data);
  data[1] = c;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) restore_solid, data);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, c);

  gts_object_destroy (GTS_OBJECT (status));
}

/**
 * gfs_cell_cm:
 * @cell: a #FttCell.
 * @cm: a #FttVector.
 *
 * Fills @cm with the coordinates of the center of mass of @cell.
 */
void gfs_cell_cm (const FttCell * cell, FttVector * cm)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (cm != NULL);

  if (GFS_IS_MIXED (cell))
    *cm = GFS_STATE (cell)->solid->cm;
  else
    ftt_cell_pos (cell, cm);
}

/**
 * gfs_solid_normal:
 * @cell: a #FttCell.
 * @n: a #FttVector.
 *
 * Fills @n with the components of the average unit normal to the
 * fraction of solid boundary contained in @cell, multiplied by the
 * area of the fraction of solid boundary contained in @cell.
 */
void gfs_solid_normal (const FttCell * cell, FttVector * n)
{
  GfsSolidVector * s;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (n != NULL);

  if ((s = GFS_STATE (cell)->solid)) {
    FttComponent c;

#if (FTT_2D)
    n->z = 0.;
#endif /* 2D */

    for (c = 0; c < FTT_DIMENSION; c++)
      (&n->x)[c] = (s->s[2*c + 1] - s->s[2*c]);
  }
  else
    n->x = n->y = n->z = 0.;
}

/**
 * gfs_face_ca:
 * @face: a #FttCellFace.
 * @ca: a #FttVector.
 *
 * Fills @ca with the coordinates of the center of area of @face.
 */
void gfs_face_ca (const FttCellFace * face, FttVector * ca)
{
  gdouble f;

  g_return_if_fail (face != NULL);
  g_return_if_fail (ca != NULL);

  ftt_face_pos (face, ca);
  if ((f = GFS_FACE_FRACTION (face)) < 1.) {
    GfsSolidVector * s = GFS_STATE (face->cell)->solid;
    gdouble h = ftt_cell_size (face->cell);
#if FTT_2D
    FttComponent cp = FTT_ORTHOGONAL_COMPONENT (face->d/2);

    (&ca->x)[cp] += (s->s[2*cp] > s->s[2*cp + 1]) ? (1. - f)/2.*h : (f - 1.)/2.*h;
#else /* 3D */
    static guint perpendicular[FTT_DIMENSION][2] = {
      {FTT_Y, FTT_Z}, {FTT_Z, FTT_X}, {FTT_X, FTT_Y}
    };
    FttComponent c0 = face->d/2;
    FttComponent c1 = perpendicular[c0][0];
    FttComponent c2 = perpendicular[c0][1];
    gboolean s1, s2;
    FttVector m, p;
    gdouble n, alpha;

    m.x = s->s[2*c1 + 1] - s->s[2*c1];
    m.y = s->s[2*c2 + 1] - s->s[2*c2];
    s1 = (m.x < 0.);
    s2 = (m.y < 0.);
    m.x = fabs (m.x);
    m.y = fabs (m.y);
    n = m.x + m.y;
    if (n > 0.) {
      m.x /= n;
      m.y /= n;
      alpha = gfs_line_alpha (&m, f);
      gfs_line_center (&m, alpha, f, &p);
      if (s1) p.x = 1. - p.x;
      if (s2) p.y = 1. - p.y;
      (&ca->x)[c1] += (p.x - 0.5)*h;
      (&ca->x)[c2] += (p.y - 0.5)*h;
    }
#endif /* 3D */
  }
}

#if !FTT_2D
static void outer_fractions_coarse_fine (FttCell * parent, FttDirection d)
{
  GfsSolidVector * solid = GFS_STATE (parent)->solid;
  FttComponent c1 = d < 4 ? 2 : 0, c2 = d < 2 || d > 3 ? 1 : 0;
  gdouble nm;
  FttVector m;
    
  m.x = solid->s[2*c1 + 1] - solid->s[2*c1]; nm = fabs (m.x);
  m.y = solid->s[2*c2 + 1] - solid->s[2*c2]; nm += fabs (m.y);
  if (nm > 0.) {
    m.x /= nm;
    m.y /= nm;
  }
  else
    m.x = 1.;
  gdouble alpha = gfs_line_alpha (&m, solid->s[d]);
  gdouble ss = 0.;
    
  FttCellChildren child;
  guint i, n = ftt_cell_children_direction (parent, d, &child);
  for (i = 0; i < n; i++)
    if (child.c[i]) {
      if (GFS_IS_MIXED (child.c[i])) {
	GfsSolidVector * s = GFS_STATE (child.c[i])->solid;
	gdouble alpha1 = alpha;
	FttVector p;
	
	ftt_cell_relative_pos (child.c[i], &p);
	alpha1 -= m.x*(0.25 + (&p.x)[c1]);
	alpha1 -= m.y*(0.25 + (&p.x)[c2]);
	
	s->s[d] = gfs_line_area (&m, 2.*alpha1);
	ss += s->s[d];
      }
      else
	ss += 1.;
    }
  /* fixme: this should not happen 
   * It happens in configurations where children cells are not cut by
   * the VOF approximation but should have non-zero surface
   * fractions */
  if (fabs (solid->s[d] - ss/n) > 1e-5)
    g_warning ("inconsistent surface fractions %d %f %f %f\n", d, solid->s[d], ss/n,
	       fabs (solid->s[d] - ss/n));
}
#endif /* 3D */

/**
 * gfs_solid_coarse_fine:
 * @parent: a mixed #FttCell with children.
 * @domain: a #GfsDomain.
 *
 * Fills the solid properties of the children of @parent.
 * Destroys all children entirely contained in the solid.
 */
void gfs_solid_coarse_fine (FttCell * parent, GfsDomain * domain)
{
  g_return_if_fail (parent);
  g_return_if_fail (domain);
  g_return_if_fail (GFS_IS_MIXED (parent));
  g_return_if_fail (!FTT_CELL_IS_LEAF (parent));

  GfsSolidVector * solid = GFS_STATE (parent)->solid;
  FttVector m;
  FttComponent c;
  gdouble n = 0;

  for (c = 0; c < FTT_DIMENSION; c++) {
    (&m.x)[c] = solid->s[2*c + 1] - solid->s[2*c];
    n += fabs ((&m.x)[c]);
  }
  if (n > 0.)
    for (c = 0; c < FTT_DIMENSION; c++)
      (&m.x)[c] /= n;
  else
    m.x = 1.;
  gdouble alpha = gfs_plane_alpha (&m, solid->a);

  gdouble h = ftt_cell_size (parent)/2.;
  guint level = ftt_cell_level (parent) + 1;
  FttCellChildren child;
  guint i;
  ftt_cell_children (parent, &child);
  for (i = 0; i < FTT_CELLS; i++) {
    gdouble alpha1 = alpha;
    FttVector p;

    ftt_cell_relative_pos (child.c[i], &p);
    for (c = 0; c < FTT_DIMENSION; c++)
      alpha1 -= (&m.x)[c]*(0.25 + (&p.x)[c]);

    if (GFS_STATE (child.c[i])->solid) {
      g_free (GFS_STATE (child.c[i])->solid);
      GFS_STATE (child.c[i])->solid = NULL;
    }

    gdouble a = gfs_plane_volume (&m, 2.*alpha1);
    if (a > 0. && a < 1.) {
      GfsSolidVector * s = GFS_STATE (child.c[i])->solid = g_malloc (sizeof (GfsSolidVector));
      s->a = a;

      ftt_cell_pos (child.c[i], &p);
      gfs_plane_center (&m, 2.*alpha1, a, &s->cm);
      gfs_plane_area_center (&m, 2.*alpha1, &s->ca);
      for (c = 0; c < FTT_DIMENSION; c++) {
	(&s->cm.x)[c] = (&p.x)[c] + h*((&s->cm.x)[c] - 0.5);
	(&s->ca.x)[c] = (&p.x)[c] + h*((&s->ca.x)[c] - 0.5);
      }

      FttDirection d;
      FttCellNeighbors n;
      ftt_cell_neighbors (child.c[i], &n);
      for (d = 0; d < FTT_NEIGHBORS; d++)
	if (!n.c[d])
	  s->s[d] = 0.;
	else if (GFS_IS_MIXED (n.c[d]) && ftt_cell_level (n.c[d]) == level)
	  s->s[d] = GFS_STATE (n.c[d])->solid->s[FTT_OPPOSITE_DIRECTION (d)];
	else if (!ftt_cell_neighbor_is_brother (child.c[i], d) && GFS_IS_FLUID (n.c[d]))
	  s->s[d] = 1.;
	else {
#if FTT_2D
	  gdouble f;
	  FttComponent c1 = d > 1, c2 = !c1;

	  if ((&m.x)[c2] == 0.) f = 0.;
	  else {
	    f = (2.*alpha1 - (&m.x)[c1]*!(d % 2))/(&m.x)[c2];
	    if (f < 0.) f = 0.; else if (f > 1.) f = 1.;
	    if ((&m.x)[c2] < 0.)
	      f = 1. - f;
	  }
	  s->s[d] = f;
#else /* 3D */
	  /* only initialises "inner" fractions */
	  if (ftt_cell_neighbor_is_brother (child.c[i], d)) {
	    FttComponent c1 = (d/2 + 1) % 3, c2 = (d/2 + 2) % 3;
	    FttVector mp;
	    mp.x = (&m.x)[c1]; 
	    mp.y = (&m.x)[c2];
	    s->s[d] = gfs_line_area (&mp, d % 2 ? 2.*alpha1 : 2.*alpha1 - (&m.x)[d/2]);
	  }
#endif /* 3D */
	}
    }
    else if (a == 0.)
      ftt_cell_destroy (child.c[i], (FttCellCleanupFunc) gfs_cell_cleanup, domain);
  }

#if !FTT_2D
  FttCellNeighbors neighbor;
  FttDirection d;
  ftt_cell_neighbors (parent, &neighbor);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (neighbor.c[d] && FTT_CELL_IS_LEAF (neighbor.c[d]) && !GFS_IS_FLUID (neighbor.c[d]))
      outer_fractions_coarse_fine (parent, d);
#endif /* 3D */
}

/**
 * Solid boundaries.
 * \beginobject{GfsSolid}
 */

static void gfs_solid_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_solid_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  gfs_generic_surface_read (GFS_SOLID (*o)->s, gfs_object_simulation (*o), fp);
}

static void gfs_solid_write (GtsObject * o, FILE * fp)
{
  GfsSimulation * sim = gfs_object_simulation (o);
  if (sim->output_solid) {
    (* GTS_OBJECT_CLASS (gfs_solid_class ())->parent_class->write) (o, fp);
    gfs_generic_surface_write (GFS_SOLID (o)->s, sim, fp);
  }
}

static void gfs_solid_destroy (GtsObject * object)
{
  gts_object_destroy (GTS_OBJECT (GFS_SOLID (object)->s));

  (* GTS_OBJECT_CLASS (gfs_solid_class ())->parent_class->destroy) (object);
}

static void gfs_solid_class_init (GtsObjectClass * klass)
{
  klass->read = gfs_solid_read;
  klass->write = gfs_solid_write;
  klass->destroy = gfs_solid_destroy;
}

static void gfs_solid_init (GfsSolid * object)
{
  object->s = GFS_GENERIC_SURFACE (gts_object_new (GTS_OBJECT_CLASS (gfs_surface_class ())));
  GFS_EVENT (object)->istep = G_MAXINT/2;
}

GfsEventClass * gfs_solid_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_solid_info = {
      "GfsSolid",
      sizeof (GfsSolid),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_solid_class_init,
      (GtsObjectInitFunc) gfs_solid_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_solid_info);
  }

  return klass;
}

/** \endobject{GfsSolid} */
