/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2013 National Institute of Water and Atmospheric Research
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
 * \brief Source terms.
 */

#include <stdlib.h>
#include <math.h>
#include "source.h"
#include "simulation.h"
#include "solid.h"
#include "init.h"
#include "vof.h"

/**
 * gfs_variable_mac_source:
 * @v: a #GfsVariable.
 * @cell: a #FttCell.
 *
 * Returns: the sum of all the sources for variable @v in @cell.
 */
gdouble gfs_variable_mac_source (GfsVariable * v, FttCell * cell)
{
  gdouble sum;
  GSList * i;

  g_return_val_if_fail (v != NULL, 0.);
  g_return_val_if_fail (cell != NULL, 0.);

  if (v->sources == NULL)
    return 0.;

  sum = 0.;
  i = GTS_SLIST_CONTAINER (v->sources)->items;
  while (i) {
    GfsSourceGeneric * s = i->data;

    if (s->mac_value)
      sum += (* s->mac_value) (s, cell, v);
    i = i->next;
  }
  return sum;
}

typedef struct {
  GfsVariable * v, * sv;
  gdouble dt;
} SourcePar;

static void add_sources (FttCell * cell, SourcePar * p)
{
  GSList * i = GTS_SLIST_CONTAINER (p->v->sources)->items;
  gdouble sum = 0;
  
  while (i) {
    GfsSourceGeneric * s = i->data;

    if (s->centered_value)
      sum += (* s->centered_value) (s, cell, p->v);
    i = i->next;
  }
  GFS_VALUE (cell, p->sv) += p->dt*sum;
}

/**
 * gfs_domain_variable_centered_sources:
 * @domain: a #GfsDomain.
 * @v: a #GfsVariable.
 * @sv: a #GfsVariable.
 * @dt: the timestep.
 *
 * Adds the source terms for @v to @sv.
 */
void gfs_domain_variable_centered_sources (GfsDomain * domain, 
					   GfsVariable * v,
					   GfsVariable * sv,
					   gdouble dt)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (v != NULL);
  g_return_if_fail (sv != NULL);

  if (v->sources) {
    SourcePar p;
    p.v = v;
    p.sv = sv;
    p.dt = dt;
    gfs_domain_cell_traverse (domain, 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) add_sources, &p);
  }
}

/**
 * gfs_domain_variable_fluxes:
 * @domain: a #GfsDomain.
 * @v: a #GfsVariable.
 * @dt: the timestep.
 *
 * Returns: a new temporary variable containing the fluxes or %NULL.
 */
GfsVariable * gfs_domain_variable_fluxes (GfsDomain * domain,
					  GfsVariable * v,
					  gdouble dt)
{
  GfsVariable * sv = NULL;

  g_return_val_if_fail (domain != NULL, NULL);
  g_return_val_if_fail (v != NULL, NULL);

  if (!v->sources)
    return NULL;

  GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;
  while (i) {
    if (GFS_SOURCE_GENERIC (i->data)->flux) {
      if (sv == NULL) {
	sv = gfs_temporary_variable (domain);
	gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) gfs_cell_reset, sv);
      }
      (* GFS_SOURCE_GENERIC (i->data)->flux) (i->data, domain, v, sv, dt);
    }
    i = i->next;
  }
  return sv;
}

/**
 * Abstract class for source terms.
 * \beginobject{GfsSourceGeneric}
 */

static void source_generic_init (GfsSourceGeneric * s)
{
  GFS_EVENT (s)->istep = 1;
}

GfsSourceGenericClass * gfs_source_generic_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_generic_info = {
      "GfsSourceGeneric",
      sizeof (GfsSourceGeneric),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) source_generic_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_source_generic_info);
  }

  return klass;
}

/**
 * gfs_source_find:
 * @v: a #GfsVariable.
 * @klass: a #GfsSourceGenericClass.
 *
 * Returns: the first source of @v descendant of @klass, or %NULL if
 * none was found.
 */
GfsSourceGeneric * gfs_source_find (GfsVariable * v, GfsSourceGenericClass * klass)
{
  g_return_val_if_fail (v != NULL, NULL);
  g_return_val_if_fail (klass != NULL, NULL);

  if (v->sources) {
    GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;
    
    while (i) {
      GtsObject * o = i->data;
      
      if (gts_object_is_from_class (o, klass))
        return GFS_SOURCE_GENERIC (o);
      i = i->next;
    }
  }
  return NULL;
}

/** \endobject{GfsSourceGeneric} */

/**
 * Source terms for scalar variables.
 * \beginobject{GfsSourceScalar}
 */

static void source_scalar_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_scalar_class ())->parent_class->write) (o, fp);

  GfsSourceScalar * s = GFS_SOURCE_SCALAR (o);
  if (s->v)
    fprintf (fp, " %s", s->v->name);
  else {
    gint i;
    fputs (" (", fp);
    for (i = 0; i < FTT_DIMENSION - 1; i++)
      fprintf (fp, "%s,", s->vector[i]->name);
    fprintf (fp, "%s)", s->vector[FTT_DIMENSION - 1]->name);
  }
}

static void source_scalar_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceScalar * source;
  GfsDomain * domain;

  if (GTS_OBJECT_CLASS (gfs_source_scalar_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_scalar_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  source = GFS_SOURCE_SCALAR (*o);
  domain =  GFS_DOMAIN (gfs_object_simulation (source));

  if (fp->type == '(') {
    /* vector */
    if (!gfs_read_variable_vector (fp, source->vector, domain))
      return;
    gint i;
    for (i = 0; i < FTT_DIMENSION; i++) {
      if (source->vector[i]->sources == NULL)
	source->vector[i]->sources = 
	  gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ()));
      gts_container_add (source->vector[i]->sources, GTS_CONTAINEE (source));
    }
  }
  else if (fp->type == GTS_STRING) {
    /* scalar */
    source->v = gfs_variable_from_name (domain->variables, fp->token->str);
    if (source->v == NULL) {
      gts_file_error (fp, "unknown variable `%s'", fp->token->str);
      return;
    }
    if (source->v->sources == NULL)
      source->v->sources = 
	gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ()));
    gts_container_add (source->v->sources, GTS_CONTAINEE (source));
    gts_file_next_token (fp);
  }
  else {
    gts_file_error (fp, "expecting a scalar or vector");
    return;
  }
}

static void source_scalar_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read =  source_scalar_read;
  GTS_OBJECT_CLASS (klass)->write = source_scalar_write;
}

GfsSourceGenericClass * gfs_source_scalar_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_scalar_info = {
      "GfsSourceScalar",
      sizeof (GfsSourceScalar),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_scalar_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_generic_class ()),
				  &gfs_source_scalar_info);
  }

  return klass;
}

/** \endobject{GfsSourceScalar} */

/**
 * Source terms for velocity.
 * \beginobject{GfsSourceVelocity}
 */

static void source_velocity_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceVelocity * source;
  GfsDomain * domain;
  FttComponent c;

  if (GTS_OBJECT_CLASS (gfs_source_velocity_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_velocity_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  source = GFS_SOURCE_VELOCITY (*o);
  domain =  GFS_DOMAIN (gfs_object_simulation (source));
  if (!(source->v = gfs_domain_velocity (domain))) {
    gts_file_error (fp, "cannot find velocity components");
    return;
  }
  for (c = 0; c < FTT_DIMENSION; c++) {
    if (source->v[c]->sources == NULL)
      source->v[c]->sources = 
	gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ()));
    gts_container_add (source->v[c]->sources, GTS_CONTAINEE (source));
  }
}

static void source_velocity_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read =  source_velocity_read;
}

GfsSourceGenericClass * gfs_source_velocity_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_info = {
      "GfsSourceVelocity",
      sizeof (GfsSourceVelocity),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_velocity_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_generic_class ()),
				  &source_info);
  }

  return klass;
}

/** \endobject{GfsSourceVelocity} */

/**
 * Source defined by a function.
 * \beginobject{GfsSource}
 */

static void source_destroy (GtsObject * o)
{
  GfsSource * s = GFS_SOURCE (o);
  if (s->intensity)
    gts_object_destroy (GTS_OBJECT (s->intensity));

  gint i;
  for (i = 0; i < FTT_DIMENSION; i++)
    if (s->intensity_v[i])
      gts_object_destroy (GTS_OBJECT (s->intensity_v[i]));

  (* GTS_OBJECT_CLASS (gfs_source_class ())->parent_class->destroy) (o);
}

static gdouble source_face_value (GfsSourceGeneric * s, 
				  FttCellFace * face, 
				  GfsVariable * v)
{
  return gfs_function_face_value (GFS_SOURCE (s)->intensity, face);
}

static gdouble source_vector_face_value (GfsSourceGeneric * s, 
					 FttCellFace * face, 
					 GfsVariable * v)
{
  FttVector p, u;
  int i;
  for (i = 0; i < FTT_DIMENSION; i++)
    (&u.x)[i] = gfs_function_face_value (GFS_SOURCE (s)->intensity_v[i], face);
  ftt_face_pos (face, &p);
  gfs_simulation_map_vector (GFS_SIMULATION (v->domain), &p, &u);
  return (&u.x)[v->component];
}

static gdouble source_value (GfsSourceGeneric * s, 
			     FttCell * cell, 
			     GfsVariable * v)
{
  return gfs_function_value (GFS_SOURCE (s)->intensity, cell);
}

static gdouble source_vector_value (GfsSourceGeneric * s, 
				    FttCell * cell, 
				    GfsVariable * v)
{
  return gfs_function_value (GFS_SOURCE (s)->intensity_v[v->component], cell);
}

static void source_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_source_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsSourceScalar * s = GFS_SOURCE_SCALAR (*o);
  if (s->v) {    
    /* scalar */
    GFS_SOURCE (s)->intensity = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_set_units (GFS_SOURCE (s)->intensity, s->v->units);
    gfs_function_read (GFS_SOURCE (s)->intensity, gfs_object_simulation (s), fp);
    if (fp->type != GTS_ERROR) {
      GfsSourceGeneric * s = GFS_SOURCE_GENERIC (*o);
      if (GFS_SOURCE_SCALAR (s)->v->face_source) {
	s->mac_value = s->centered_value = NULL;
	s->face_value = source_face_value;
      }
    }
  }
  else {
    /* vector */
    if (!gfs_read_function_vector (fp, s->vector, GFS_SOURCE (s)->intensity_v, 
				   gfs_object_simulation (s)))
      return;
    if (s->vector[0]->face_source) {
      GfsSourceGeneric * s = GFS_SOURCE_GENERIC (*o);
      s->mac_value = s->centered_value = NULL;
      s->face_value = source_vector_face_value;
    }
    else {
      GfsSourceGeneric * s = GFS_SOURCE_GENERIC (*o);
      s->mac_value = s->centered_value = source_vector_value;
    }
    gts_file_next_token (fp);
  }
}

static void source_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_source_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_source_class ())->parent_class->write) 
      (o, fp);
  GfsSource * s = GFS_SOURCE (o);
  if (s->intensity)
    gfs_function_write (s->intensity, fp);
  else {
    gint i;
    fputs (" (", fp);
    for (i = 0; i < FTT_DIMENSION - 1; i++) {
      gfs_function_write (s->intensity_v[i], fp);
      fputc (',', fp);
    }
    gfs_function_write (s->intensity_v[FTT_DIMENSION - 1], fp);
    fputc (')', fp);
  }    
}

static void source_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = source_destroy;
  GTS_OBJECT_CLASS (klass)->read = source_read;
  GTS_OBJECT_CLASS (klass)->write = source_write;
}

static void source_init (GfsSourceGeneric * s)
{
  s->mac_value = s->centered_value = source_value;
}

GfsSourceGenericClass * gfs_source_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_info = {
      "GfsSource",
      sizeof (GfsSource),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_class_init,
      (GtsObjectInitFunc) source_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_scalar_class ()),
				  &source_info);
  }

  return klass;
}

/** \endobject{GfsSource} */

/**
 * Controlling the spatially-averaged value of a scalar.
 * \beginobject{GfsSourceControl}
 */

static void source_control_destroy (GtsObject * o)
{
  if (GFS_SOURCE_CONTROL (o)->intensity)
    gts_object_destroy (GTS_OBJECT (GFS_SOURCE_CONTROL (o)->intensity));

  (* GTS_OBJECT_CLASS (gfs_source_control_class ())->parent_class->destroy) (o);
}

static void source_control_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_control_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GFS_SOURCE_CONTROL (*o)->intensity = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_set_units (GFS_SOURCE_CONTROL (*o)->intensity, GFS_SOURCE_SCALAR (*o)->v->units);
  gfs_function_read (GFS_SOURCE_CONTROL (*o)->intensity, gfs_object_simulation (*o), fp);
}

static void source_control_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_control_class ())->parent_class->write) (o, fp);
  gfs_function_write (GFS_SOURCE_CONTROL (o)->intensity, fp);
}

typedef struct {
  GfsVariable * v;
  gdouble s, sv;
} Sum;

static void sum (FttCell * cell, Sum * s)
{
  gdouble vol = gfs_cell_volume (cell, s->v->domain);
  s->s += vol*GFS_VALUE (cell, s->v);
  s->sv += vol;
}

static gboolean source_control_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* gfs_event_class ()->event) (event, sim)) {
    GfsSourceControl * s = GFS_SOURCE_CONTROL (event);
    GfsDomain * domain = GFS_DOMAIN (sim);
    Sum su = { GFS_SOURCE_SCALAR (event)->v, 0., 0. };
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) sum, &su);
    gfs_all_reduce (domain, su.s, MPI_DOUBLE, MPI_SUM);
    gfs_all_reduce (domain, su.sv, MPI_DOUBLE, MPI_SUM);
    s->s = sim->advection_params.dt > 0. && su.sv > 0. ? 
      (gfs_function_value (s->intensity, NULL) - su.s/su.sv)/sim->advection_params.dt: 0.;
    return TRUE;
  }
  return FALSE;
}

static void source_control_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = source_control_read;
  GTS_OBJECT_CLASS (klass)->write = source_control_write;
  GTS_OBJECT_CLASS (klass)->destroy = source_control_destroy;
  GFS_EVENT_CLASS (klass)->event = source_control_event;
}

static gdouble source_control_value (GfsSourceGeneric * s, 
				     FttCell * cell, 
				     GfsVariable * v)
{
  return GFS_SOURCE_CONTROL (s)->s;
}

static void source_control_init (GfsSourceGeneric * s)
{
  s->mac_value = s->centered_value = source_control_value;
}

GfsSourceGenericClass * gfs_source_control_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_control_info = {
      "GfsSourceControl",
      sizeof (GfsSourceControl),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_control_class_init,
      (GtsObjectInitFunc) source_control_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_scalar_class ()),
				  &source_control_info);
  }

  return klass;
}

/** \endobject{GfsSourceControl} */

/**
 * Controlling the value of a scalar averaged at a given spatial scale.
 * \beginobject{GfsSourceControlField}
 */

static void source_control_field_destroy (GtsObject * o)
{
  if (GFS_SOURCE_CONTROL_FIELD (o)->s)
    gts_object_destroy (GTS_OBJECT (GFS_SOURCE_CONTROL_FIELD (o)->s));

  (* GTS_OBJECT_CLASS (gfs_source_control_field_class ())->parent_class->destroy) (o);
}

static gdouble source_control_field_face_value (GfsSourceGeneric * s, 
						FttCellFace * face, 
						GfsVariable * v)
{
  return gfs_face_interpolated_value_generic (face, GFS_SOURCE_CONTROL_FIELD (s)->s);
}

static void source_control_field_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_control_field_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsSourceControlField * f = GFS_SOURCE_CONTROL_FIELD (*o);
  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting an integer (level)");
    return;
  }
  f->level = atoi (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_DOUBLE) {
    gts_file_error (fp, "expecting a number (tau)");
    return;
  }
  f->tau = atof (fp->token->str);
  if (f->tau <= 0.) {
    gts_file_error (fp, "tau must be strictly positive");
    return;
  }
  gts_file_next_token (fp);  

  f->s = gfs_temporary_variable (GFS_DOMAIN (gfs_object_simulation (*o)));

  GfsSourceGeneric * s = GFS_SOURCE_GENERIC (*o);
  if (GFS_SOURCE_SCALAR (s)->v->face_source) {
    s->mac_value = s->centered_value = NULL;
    s->face_value = source_control_field_face_value;
  }
}

static void source_control_field_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_control_field_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %d %g", GFS_SOURCE_CONTROL_FIELD (o)->level, GFS_SOURCE_CONTROL_FIELD (o)->tau);
}

static void source_control_field_root (FttCell * root, GfsSourceControlField * f)
{
  Sum su = { GFS_SOURCE_SCALAR (f)->v, 0., 0. };
  ftt_cell_traverse (root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		     (FttCellTraverseFunc) sum, &su);
  GFS_VALUE (root, f->s) = su.sv > 0. ? 
    (gfs_function_value (GFS_SOURCE_CONTROL (f)->intensity, root) - su.s/su.sv)/f->tau: 0.;
}

typedef struct {
  FttCell * root;
  gdouble * corners;
  GfsVariable * v;
} ExtraData;

static void extrapolate (FttCell * cell, ExtraData * p)
{
  FttVector pos;
  ftt_cell_pos (cell, &pos);
  GFS_VALUE (cell, p->v) = gfs_interpolate_from_corners (p->root, pos, p->corners);
}

static void extrapolate_field (FttCell * root, GfsSourceControlField * f)
{
  gdouble corners[4*(FTT_DIMENSION - 1) + 1];
  ExtraData p = { root, corners, f->s };
  gfs_cell_corner_values (root, f->s, f->level, corners);
  g_assert (!FTT_CELL_IS_LEAF (root));
  ftt_cell_traverse (root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		     (FttCellTraverseFunc) extrapolate, &p);
}

static gboolean source_control_field_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* gfs_event_class ()->event) (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsSourceControlField * f = GFS_SOURCE_CONTROL_FIELD (event);
    gfs_catch_floating_point_exceptions ();
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, 
			      FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS, f->level,
			      (FttCellTraverseFunc) source_control_field_root, event);
    gfs_restore_fpe_for_function (GFS_SOURCE_CONTROL (event)->intensity);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS, f->level, f->s);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER,
			      FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_NON_LEAFS, f->level,
			      (FttCellTraverseFunc) extrapolate_field, event);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, f->s);
    return TRUE;
  }
  return FALSE;
}

static void source_control_field_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = source_control_field_destroy;
  GTS_OBJECT_CLASS (klass)->read = source_control_field_read;
  GTS_OBJECT_CLASS (klass)->write = source_control_field_write;
  GFS_EVENT_CLASS (klass)->event = source_control_field_event;
}

static gdouble source_control_field_value (GfsSourceGeneric * s, 
					   FttCell * cell, 
					   GfsVariable * v)
{
  return GFS_VALUE (cell, GFS_SOURCE_CONTROL_FIELD (s)->s);
}

static void source_control_field_init (GfsSourceGeneric * s)
{
  s->mac_value = NULL;
  s->centered_value = source_control_field_value;
  GFS_SOURCE_CONTROL_FIELD (s)->tau = 1.;
}

GfsSourceGenericClass * gfs_source_control_field_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_control_field_info = {
      "GfsSourceControlField",
      sizeof (GfsSourceControlField),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_control_field_class_init,
      (GtsObjectInitFunc) source_control_field_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_control_class ()),
				  &source_control_field_info);
  }

  return klass;
}

/** \endobject{GfsSourceControlField} */

/**
 * Adding a flux on a given area or volume.
 * \beginobject{GfsSourceFlux}
 */

static void source_flux_destroy (GtsObject * o)
{
  if (GFS_SOURCE_FLUX (o)->intensity)
    gts_object_destroy (GTS_OBJECT (GFS_SOURCE_FLUX (o)->intensity));
  if (GFS_SOURCE_FLUX (o)->fraction)
    gts_object_destroy (GTS_OBJECT (GFS_SOURCE_FLUX (o)->fraction));

  (* GTS_OBJECT_CLASS (gfs_source_flux_class ())->parent_class->destroy) (o);
}

static void source_flux_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_flux_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsSourceFlux * f = GFS_SOURCE_FLUX (*o);
  f->intensity = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_set_units (f->intensity, GFS_SOURCE_SCALAR (f)->v->units + FTT_DIMENSION);
  gfs_function_read (f->intensity, gfs_object_simulation (f), fp);
  if (fp->type == GTS_ERROR)
    return;
  f->fraction = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_read (f->fraction, gfs_object_simulation (f), fp);
}

static void source_flux_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_flux_class ())->parent_class->write) (o, fp);
  gfs_function_write (GFS_SOURCE_FLUX (o)->intensity, fp);
  gfs_function_write (GFS_SOURCE_FLUX (o)->fraction, fp);
}

static void add (FttCell * cell, GfsSourceFlux * s)
{
  s->s += gfs_cell_volume (cell, GFS_DOMAIN (gfs_object_simulation (s)))*
    gfs_function_value (s->fraction, cell);
}

static gboolean source_flux_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* gfs_event_class ()->event) (event, sim)) {
    GfsSourceFlux * s = GFS_SOURCE_FLUX (event);
    s->s = 0.;
    gfs_catch_floating_point_exceptions ();
    gfs_domain_traverse_leaves (GFS_DOMAIN (sim), (FttCellTraverseFunc) add, s);
    gfs_restore_fpe_for_function (s->fraction);
    gfs_catch_floating_point_exceptions ();
    s->s = s->s > 0. ? gfs_function_value (s->intensity, NULL)/s->s : 0.;
    gfs_restore_fpe_for_function (s->intensity);
    return TRUE;
  }
  return FALSE;
}

static void source_flux_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = source_flux_read;
  GTS_OBJECT_CLASS (klass)->write = source_flux_write;
  GTS_OBJECT_CLASS (klass)->destroy = source_flux_destroy;
  GFS_EVENT_CLASS (klass)->event = source_flux_event;
}

static gdouble source_flux_value (GfsSourceGeneric * s, 
				  FttCell * cell, 
				  GfsVariable * v)
{
  return GFS_SOURCE_FLUX (s)->s*gfs_function_value (GFS_SOURCE_FLUX (s)->fraction, cell);
}

static void source_flux_init (GfsSourceGeneric * s)
{
  s->mac_value = s->centered_value = source_flux_value;
}

GfsSourceGenericClass * gfs_source_flux_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_flux_info = {
      "GfsSourceFlux",
      sizeof (GfsSourceFlux),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_flux_class_init,
      (GtsObjectInitFunc) source_flux_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_scalar_class ()),
				  &source_flux_info);
  }

  return klass;
}

/** \endobject{SourceFlux} */

/**
 * Diffusion coefficient.
 * \beginobject{GfsDiffusion}
 */

static void diffusion_destroy (GtsObject * o)
{
  GfsDiffusion * d = GFS_DIFFUSION (o);

  if (d->mu && d->mu != gfs_function_get_variable (d->val))
    gts_object_destroy (GTS_OBJECT (d->mu));
  gts_object_destroy (GTS_OBJECT (d->val));

  (* GTS_OBJECT_CLASS (gfs_diffusion_class ())->parent_class->destroy) (o);
}

static void diffusion_read (GtsObject ** o, GtsFile * fp)
{
  GfsDiffusion * d = GFS_DIFFUSION (*o);

  gfs_function_read (d->val, gfs_object_simulation (*o), fp);
  if (fp->type == GTS_ERROR)
    return;
  gfs_function_set_units (d->val, 2.);

  if (fp->type == '{') {
    gfs_multilevel_params_read (&d->par, fp);
    if (fp->type == GTS_ERROR)
    return;
  }
}

static void diffusion_write (GtsObject * o, FILE * fp)
{
  gfs_function_write (GFS_DIFFUSION (o)->val, fp);
  fputc (' ', fp);
  gfs_multilevel_params_write (&GFS_DIFFUSION (o)->par, fp);
}

static void update_mu (FttCell * cell, GfsDiffusion * d)
{
  GFS_VALUE (cell, d->mu) = gfs_function_value (d->val, cell);
}

static gboolean diffusion_event (GfsEvent * event, GfsSimulation * sim)
{
  GfsDiffusion * d = GFS_DIFFUSION (event);

  if (!gfs_function_is_constant (d->val)) {
    if (d->mu == NULL && (d->mu = gfs_function_get_variable (d->val)) == NULL)
      d->mu = gfs_domain_add_variable (GFS_DOMAIN (sim), NULL, NULL);
    if (d->mu != gfs_function_get_variable (d->val)) {
      gfs_catch_floating_point_exceptions ();
      gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) update_mu, event);
      gfs_restore_fpe_for_function (d->val);
    }
    gfs_domain_cell_traverse (GFS_DOMAIN (sim),
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_get_from_below_intensive, d->mu);
    gfs_domain_bc (GFS_DOMAIN (sim), FTT_TRAVERSE_ALL, -1, d->mu);
    return TRUE;
  }
  return FALSE;
}

static gdouble diffusion_face (GfsDiffusion * d, FttCellFace * face)
{
  if (d->par.function) return gfs_function_face_value (d->val, face);
  if (d->mu) return gfs_face_interpolated_value_generic (face, d->mu);
  gdouble val = gfs_function_get_constant_value (d->val);
  return val < G_MAXDOUBLE ? val : 0.;
}

static gdouble diffusion_cell (GfsDiffusion * d, FttCell * cell)
{
  if (d->mu) return GFS_VALUE (cell, d->mu);
  gdouble val = gfs_function_get_constant_value (d->val);
  return val < G_MAXDOUBLE ? val : 0.;
}

static gdouble diffusion_concentration_face (GfsDiffusion * d, FttCellFace * face)
{
  return diffusion_face (d, face)*gfs_vof_face_fraction (face, GFS_VARIABLE_TRACER_VOF (d->phase));
}

static gdouble diffusion_concentration_cell (GfsDiffusion * d, FttCell * cell)
{
  return diffusion_cell (d, cell)*GFS_VALUE (cell, d->phase);
}

static void diffusion_class_init (GfsDiffusionClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = diffusion_destroy;
  GTS_OBJECT_CLASS (klass)->read = diffusion_read;
  GTS_OBJECT_CLASS (klass)->write = diffusion_write;
  GFS_EVENT_CLASS (klass)->event = diffusion_event;
}

static void diffusion_init (GfsDiffusion * d)
{
  gfs_multilevel_params_init (&d->par);
  d->par.tolerance = 1e-6;
  d->val = gfs_function_new (gfs_function_class (), 0.);
  d->mu = NULL;
  d->face = diffusion_face;
  d->cell = diffusion_cell;
}

GfsDiffusionClass * gfs_diffusion_class (void)
{
  static GfsDiffusionClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo diffusion_info = {
      "GfsDiffusion",
      sizeof (GfsDiffusion),
      sizeof (GfsDiffusionClass),
      (GtsObjectClassInitFunc) diffusion_class_init,
      (GtsObjectInitFunc) diffusion_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &diffusion_info);
  }

  return klass;
}

gdouble gfs_diffusion_face (GfsDiffusion * d, FttCellFace * f)
{
  return (* d->face) (d, f);
}

gdouble gfs_diffusion_cell (GfsDiffusion * d, FttCell * cell)
{
  return (* d->cell) (d, cell);
}

/** \endobject{GfsDiffusion} */

/**
 * Scalar diffusion.
 * \beginobject{GfsSourceDiffusion}
 */

static void source_diffusion_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_SOURCE_DIFFUSION (o)->D));

  (* GTS_OBJECT_CLASS (gfs_source_diffusion_class ())->parent_class->destroy) (o);
}

static GfsSourceDiffusion * previous_diffusion_source (GfsVariable * v,
						       GfsSourceDiffusion * d)
{
  GSList * i;
  i = GTS_SLIST_CONTAINER (v->sources)->items;
  while (i) {
    if (i->data != d &&
	GFS_IS_SOURCE_DIFFUSION (i->data) && 
	GFS_SOURCE_DIFFUSION (i->data)->phi == d->phi)
      return i->data;
    i = i->next;
  }
  return NULL;
}

static void source_diffusion_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceDiffusion * d;

  if (GTS_OBJECT_CLASS (gfs_source_diffusion_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_diffusion_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  d = GFS_SOURCE_DIFFUSION (*o);
  if (previous_diffusion_source (GFS_SOURCE_SCALAR (d)->v, d)) {
    gts_file_error (fp, "only one diffusion source can be specified");
    return;
  }

  gfs_object_simulation_set (d->D, gfs_object_simulation (d));
  (* GTS_OBJECT (d->D)->klass->read) ((GtsObject **) &d->D, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (fp->type != GTS_STRING)
    d->phi = GFS_SOURCE_SCALAR (*o)->v;
  else {
    d->phi = gfs_variable_from_name (domain->variables, fp->token->str);
    if (!d->phi) {
      gts_file_error (fp, "unknown variable '%s'", fp->token->str);
      return;
    }
    gts_file_next_token (fp);
    gfs_function_set_units (d->D->val, 
			    2. + GFS_SOURCE_SCALAR (*o)->v->units - d->phi->units);
  }
  if (GFS_IS_VARIABLE_VOF_CONCENTRATION (GFS_SOURCE_SCALAR (d)->v)) {
    d->D->phase = GFS_VARIABLE (GFS_VARIABLE_VOF_CONCENTRATION (GFS_SOURCE_SCALAR (d)->v)->vof);
    d->D->face = diffusion_concentration_face;
    d->D->cell = diffusion_concentration_cell;
  }

  if (GFS_SOURCE_SCALAR (d)->v->component < FTT_DIMENSION && !gfs_function_is_constant (d->D->val))
    g_warning ("%d:%d: Terms may be missing when using variable diffusion\n"
	       "on vector quantities",
	       fp->line, fp->pos);
}

static void source_diffusion_write (GtsObject * o, FILE * fp)
{
  GfsSourceDiffusion * d = GFS_SOURCE_DIFFUSION (o);

  (* GTS_OBJECT_CLASS (gfs_source_diffusion_class ())->parent_class->write) (o, fp);
  (* GTS_OBJECT (d->D)->klass->write) (GTS_OBJECT (d->D), fp);

  if (d->phi != GFS_SOURCE_SCALAR (o)->v)
    fprintf (fp, " %s", d->phi->name);
}

static gboolean source_diffusion_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* gfs_event_class ()->event) (event, sim)) {
    GfsSourceDiffusion * d = GFS_SOURCE_DIFFUSION (event);

    if ((* GFS_EVENT_CLASS (GTS_OBJECT (d->D)->klass)->event))
      (* GFS_EVENT_CLASS (GTS_OBJECT (d->D)->klass)->event) (GFS_EVENT (d->D), sim);
    return TRUE;
  }
  return FALSE;
}

static gdouble source_diffusion_value (GfsSourceGeneric * s, 
				       FttCell * cell,
				       GfsVariable * v)
{
  GfsVariable * phi = GFS_SOURCE_DIFFUSION (s)->phi != GFS_SOURCE_SCALAR (s)->v ? 
    GFS_SOURCE_DIFFUSION (s)->phi : v;
  FttCellFace f;
  FttCellNeighbors n;
  GfsGradient g = { 0., 0. };
  gdouble v0, h;

  if (GFS_IS_MIXED (cell)) /* this improves results for channel test */
    return 0.;

  v0 = GFS_VALUE (cell, phi);
  f.cell = cell;
  ftt_cell_neighbors (cell, &n);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    gdouble D;

    f.neighbor = n.c[f.d];
    D = gfs_source_diffusion_face (GFS_SOURCE_DIFFUSION (s), &f);
    if (f.neighbor) {
      GfsGradient e;

      gfs_face_gradient (&f, &e, phi->i, -1);
      g.a += D*e.a;
      g.b += D*e.b;
    }
    else if (f.d/2 == phi->component) {
      g.a += D;
      g.b -= D*v0;
    }
  }
  h = ftt_cell_size (cell);

  GfsFunction * alpha = phi->component < FTT_DIMENSION ? 
    gfs_object_simulation (s)->physical_params.alpha : NULL;
  return (alpha ? gfs_function_value (alpha, cell) : 1.)*(g.b - g.a*v0)/(h*h);
}

static void source_diffusion_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = source_diffusion_destroy;
  GTS_OBJECT_CLASS (klass)->read = source_diffusion_read;
  GTS_OBJECT_CLASS (klass)->write = source_diffusion_write;

  GFS_EVENT_CLASS (klass)->event = source_diffusion_event;
}

static void source_diffusion_init (GfsSourceDiffusion * d)
{
  d->D = GFS_DIFFUSION (gts_object_new (GTS_OBJECT_CLASS (gfs_diffusion_class ())));
  GFS_SOURCE_GENERIC (d)->mac_value = source_diffusion_value;
}

GfsSourceGenericClass * gfs_source_diffusion_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_diffusion_info = {
      "GfsSourceDiffusion",
      sizeof (GfsSourceDiffusion),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_diffusion_class_init,
      (GtsObjectInitFunc) source_diffusion_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_scalar_class ()),
				  &source_diffusion_info);
  }

  return klass;
}

gdouble gfs_source_diffusion_face (GfsSourceDiffusion * d, FttCellFace * f)
{
  g_return_val_if_fail (d != NULL, 0.);
  g_return_val_if_fail (f != NULL, 0.);

  return gfs_diffusion_face (d->D, f);
}

gdouble gfs_source_diffusion_cell (GfsSourceDiffusion * d, FttCell * cell)
{
  g_return_val_if_fail (d != NULL, 0.);
  g_return_val_if_fail (cell != NULL, 0.);

  return gfs_diffusion_cell (d->D, cell);
}

/** \endobject{GfsSourceDiffusionExplicit} */

/**
 * Time-explicit scalar diffusion.
 * \beginobject{GfsSourceDiffusionExplicit}
 */

typedef struct {
  GfsSourceGeneric * s;
  GfsVariable * v, * sv, * tv;
  gdouble dt;
} FluxPar;

static void add_diffusion_explicit_flux (FttCell * cell, FluxPar * p)
{
  /* see also poisson.c:diffusion_rhs() */
  gdouble f, h, val;
  FttCellNeighbors neighbor;
  FttCellFace face;
  
  if (GFS_IS_MIXED (cell)) {
    GfsSolidVector * solid = GFS_STATE (cell)->solid;
    if (((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
      f = gfs_cell_dirichlet_gradient_flux (cell, p->v->i, -1, solid->fv);
    else
      f = solid->fv*solid->v.x;
  }
  else
    f = 0.; /* Neumann condition by default */
  h = ftt_cell_size (cell);
  val = GFS_VALUE (cell, p->v);
  face.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++) {
    GfsGradient g;

    face.neighbor = neighbor.c[face.d];
    gfs_face_cm_weighted_gradient (&face, &g, p->v->i, -1);
    if (face.d/2 == p->v->component) {
      g.a *= 2.;
      g.b *= 2.;
    }
    f += g.b - g.a*val;
  }
  GFS_VALUE (cell, p->sv) += f/(h*h);
}

static void source_diffusion_explicit_flux (GfsSourceGeneric * s, 
					    GfsDomain * domain, 
					    GfsVariable * v, GfsVariable * sv, 
					    gdouble dt)
{
  GfsVariable * phi = GFS_SOURCE_DIFFUSION (s)->phi;
  gfs_domain_surface_bc (domain, phi);
  gfs_diffusion_coefficients (domain, GFS_SOURCE_DIFFUSION (s), dt, NULL, NULL, NULL, 1.);
  FluxPar p = { s, phi, sv };
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) add_diffusion_explicit_flux, &p);
}

static void gfs_source_diffusion_explicit_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_source_diffusion_explicit_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_diffusion_explicit_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
}

static void gfs_source_diffusion_explicit_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_diffusion_explicit_class ())->parent_class->write) (o, fp);
}

typedef struct {
  GfsFunction * alpha;
  GfsSourceGeneric * s;
  gdouble dtmax;
} StabilityParams;

static void cell_diffusion_stability (FttCell * cell,
				      StabilityParams * par)
{
  if (GFS_IS_MIXED (cell))
    return;

  FttCellFace f;
  FttCellNeighbors n;
  gdouble Dmax = 0.;
  f.cell = cell;
  ftt_cell_neighbors (cell, &n);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    gdouble D;

    f.neighbor = n.c[f.d];
    D = gfs_source_diffusion_face (GFS_SOURCE_DIFFUSION (par->s), &f);
    if (D > Dmax)
      Dmax = D;
  }

  gdouble h = ftt_cell_size (cell);
  if (Dmax > 0.) {
    gdouble dtmax = h*h/(Dmax*(par->alpha ? gfs_function_value (par->alpha, cell) : 1.))
      /* fixme: why 10? */
      /10.;
    if (dtmax < par->dtmax)
      par->dtmax = dtmax;
  }
}

static gdouble source_diffusion_stability (GfsSourceGeneric * s,
					   GfsSimulation * sim)
{
  StabilityParams par;

  par.s = s;
  par.dtmax = G_MAXDOUBLE;
  par.alpha = NULL;
  if (GFS_SOURCE_SCALAR (s)->v == GFS_SOURCE_DIFFUSION (s)->phi)
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) cell_diffusion_stability, &par);
  /* else
     fixme: what is the stability condition in that case? */
  return par.dtmax;
}

static void gfs_source_diffusion_explicit_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_source_diffusion_explicit_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_source_diffusion_explicit_write;
  klass->stability = source_diffusion_stability;
}

static void gfs_source_diffusion_explicit_init (GfsSourceGeneric * s)
{
  s->mac_value = s->centered_value = NULL;
  s->flux = source_diffusion_explicit_flux;
}

GfsSourceGenericClass * gfs_source_diffusion_explicit_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_diffusion_explicit_info = {
      "GfsSourceDiffusionExplicit",
      sizeof (GfsSourceDiffusionExplicit),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_diffusion_explicit_class_init,
      (GtsObjectInitFunc) gfs_source_diffusion_explicit_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_diffusion_class ()),
				  &gfs_source_diffusion_explicit_info);
  }

  return klass;
}

/** \endobject{GfsSourceDiffusionExplicit} */

/**
 * Viscous terms.
 * \beginobject{GfsSourceViscosity}
 */

static void compute_transverse (FttCell * cell, FluxPar * p)
{
  GfsVariable ** v = GFS_SOURCE_VISCOSITY (p->s)->v;
  gdouble h = ftt_cell_size (cell);
  GFS_VALUE (cell, p->tv) = gfs_center_gradient (cell, 
						 p->v->component, 
						 v[(p->v->component + 1) % FTT_DIMENSION]->i)/h;
}

static void source_viscosity_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceViscosity * source;
  GfsSourceDiffusion * d;
  GfsDomain * domain;
  FttComponent c;

  (* GTS_OBJECT_CLASS (gfs_source_velocity_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  source = GFS_SOURCE_VISCOSITY (*o);
  domain =  GFS_DOMAIN (gfs_object_simulation (source));
  if (!(source->v = gfs_domain_velocity (domain))) {
    gts_file_error (fp, "cannot find velocity components");
    return;
  }
  for (c = 0; c < FTT_DIMENSION; c++) {
    if (source->v[c]->sources == NULL)
      source->v[c]->sources = 
	gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ()));
    gts_container_add (source->v[c]->sources, GTS_CONTAINEE (source));
  }

  d = GFS_SOURCE_DIFFUSION (*o);
  gfs_object_simulation_set (d->D, gfs_object_simulation (d));
  (* GTS_OBJECT (d->D)->klass->read) ((GtsObject **) &d->D, fp);
}

static void source_viscosity_write (GtsObject * o, FILE * fp)
{
  GfsSourceDiffusion * d = GFS_SOURCE_DIFFUSION (o);

  (* GTS_OBJECT_CLASS (gfs_source_velocity_class ())->parent_class->write) (o, fp);
  (* GTS_OBJECT (d->D)->klass->write) (GTS_OBJECT (d->D), fp);
}

/* Returns explicit terms not included in the implicit scheme
   i.e. terms linked to a variable viscosity and/or metric terms */
static gdouble source_viscosity_non_diffusion_value (GfsSourceGeneric * s,
						     FttCell * cell,
						     GfsVariable * v)
{
  gdouble src = 0.;
  GfsVariable * mu = GFS_SOURCE_DIFFUSION (s)->D->mu;

  /* variable viscosity */
  if (mu != NULL) {
    GfsVariable ** u = GFS_SOURCE_VISCOSITY (s)->v;
    FttComponent c = v->component, j;
    GfsFunction * alpha = gfs_object_simulation (s)->physical_params.alpha;
    gdouble h = ftt_cell_size (cell);
    gdouble a = 0.;

    for (j = 0; j < FTT_DIMENSION; j++)
      a += (gfs_center_gradient (cell, c, u[j]->i)*
	    gfs_center_gradient (cell, j, mu->i));
    src += a*(alpha ? gfs_function_value (alpha, cell) : 1.)/(h*h);
  }

  /* metric */
  if (v->domain->viscous_metric_explicit)
    src += (* v->domain->viscous_metric_explicit) (v->domain, cell, v, GFS_SOURCE_DIFFUSION (s)->D);

  return src;
}

static gdouble source_viscosity_value (GfsSourceGeneric * s,
				       FttCell * cell,
				       GfsVariable * v)
{
  return (source_diffusion_value (s, cell, v) +
	  source_viscosity_non_diffusion_value (s, cell, v));
}

static void source_viscosity_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = source_viscosity_read;
  GTS_OBJECT_CLASS (klass)->write = source_viscosity_write;
}

static void source_viscosity_init (GfsSourceGeneric * s)
{
  s->mac_value = source_viscosity_value;
  s->centered_value = source_viscosity_non_diffusion_value;
  s->flux = NULL;
}

GfsSourceGenericClass * gfs_source_viscosity_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_viscosity_info = {
      "GfsSourceViscosity",
      sizeof (GfsSourceViscosity),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_viscosity_class_init,
      (GtsObjectInitFunc) source_viscosity_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_diffusion_class ()),
				  &source_viscosity_info);
  }

  return klass;
}

/** \endobject{GfsSourceViscosity} */

/**
 *
 * \beginobject{GfsSourceViscosityExplicit}
 */

static gdouble source_viscosity_stability (GfsSourceGeneric * s,
					   GfsSimulation * sim)
{
  StabilityParams par;

  par.s = s;
  par.dtmax = G_MAXDOUBLE;
  par.alpha = sim->physical_params.alpha;
  gfs_catch_floating_point_exceptions ();
  gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) cell_diffusion_stability, &par);
  gfs_restore_fpe_for_function (par.alpha);
  return par.dtmax;
}

static void source_viscosity_explicit_class_init (GfsSourceGenericClass * klass)
{
  klass->stability = source_viscosity_stability;
}

/* Defining this will use the divergence-free condition to decouple
   the diffusion equations for each component. This only works of
   course for constant viscosity and does not work for axisymmetric flows. */
/* #define NOTRANSVERSE */

static void add_viscosity_explicit_flux (FttCell * cell, FluxPar * p)
{
  FttCellFace f;
  FttCellNeighbors n;
  GfsGradient g = { 0., 0. };
  gdouble v0;

  if (GFS_IS_MIXED (cell)) {
    if (((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
      g.b = gfs_cell_dirichlet_gradient_flux (cell, p->v->i, -1., 0.);
  }

  v0 = GFS_VALUE (cell, p->v);
  f.cell = cell;
  ftt_cell_neighbors (cell, &n);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    GfsGradient e;

    f.neighbor = n.c[f.d];
    gfs_face_cm_weighted_gradient (&f, &e, p->v->i, -1);
#ifndef NOTRANSVERSE
    if (f.d/2 == p->v->component) {
      e.a *= 2.;
      e.b *= 2.;
    }
#endif
    g.a += e.a;
    g.b += e.b;
  }

  gdouble transverse = 0.;
#ifndef NOTRANSVERSE
#if FTT_2D
  FttComponent ortho = (p->v->component + 1) % FTT_DIMENSION;

  for (f.d = 2*ortho; f.d <= 2*ortho + 1; f.d++) {
    f.neighbor = n.c[f.d];
    transverse += (FTT_FACE_DIRECT (&f) ? 1. : -1.)*
      gfs_face_weighted_interpolated_value (&f, p->tv->i);
  }
#else
  g_assert_not_implemented ();
#endif  
#endif

  GfsFunction * alpha = gfs_object_simulation (p->s)->physical_params.alpha;
  gdouble h = ftt_cell_size (cell);
  GFS_VALUE (cell, p->sv) += (alpha ? gfs_function_value (alpha, cell) : 1.)*
    ((g.b - g.a*v0)/h + transverse)/h;
}

static void add_axisymmetric_term (FttCell * cell, FluxPar * p)
{
  GfsFunction * alpha = gfs_object_simulation (p->s)->physical_params.alpha;
  gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
  GFS_VALUE (cell, p->sv) -= 
    (alpha ? gfs_function_value (alpha, cell) : 1.)*
    2.*gfs_source_diffusion_cell (GFS_SOURCE_DIFFUSION (p->s), cell)*
    GFS_VALUE (cell, p->v)*
    a*a/gfs_domain_cell_fraction (p->v->domain, cell)*
    p->dt;
}

static void source_viscosity_explicit_flux (GfsSourceGeneric * s, 
					    GfsDomain * domain, 
					    GfsVariable * v, GfsVariable * sv, 
					    gdouble dt)
{
  FluxPar p;

  gfs_domain_surface_bc (domain, v);
  gfs_diffusion_coefficients (domain, GFS_SOURCE_DIFFUSION (s), dt, NULL, NULL, NULL, 1.);
  p.s = s;
  p.v = v;
  p.sv = sv;
  p.dt = dt;
  p.tv = gfs_temporary_variable (domain);
  gfs_traverse_and_bc (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) compute_transverse, &p,
		       p.tv, p.tv);
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) add_viscosity_explicit_flux, &p);
  if (GFS_IS_AXI (domain) && v->component == FTT_Y)
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) add_axisymmetric_term, &p);
  gts_object_destroy (GTS_OBJECT (p.tv));
}

static void source_viscosity_explicit_init (GfsSourceGeneric * s)
{
  s->mac_value = NULL;
  s->centered_value = NULL;
  s->flux = source_viscosity_explicit_flux;
}

GfsSourceGenericClass * gfs_source_viscosity_explicit_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_viscosity_explicit_info = {
      "GfsSourceViscosityExplicit",
      sizeof (GfsSourceViscosity),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_viscosity_explicit_class_init,
      (GtsObjectInitFunc) source_viscosity_explicit_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_viscosity_class ()),
				  &source_viscosity_explicit_info);
  }

  return klass;
}

/** \endobject{GfsSourceViscosityExplicit} */

/**
 * Coriolis acceleration, linear friction and/or advection metric terms.
 * \beginobject{GfsSourceCoriolis}
 */

static void source_coriolis_destroy (GtsObject * o)
{
  GfsSourceCoriolis * s = GFS_SOURCE_CORIOLIS (o);
  FttComponent c;

  if (s->omegaz)
    gts_object_destroy (GTS_OBJECT (s->omegaz));
  if (s->drag)
    gts_object_destroy (GTS_OBJECT (s->drag));

  for (c = 0; c <  FTT_DIMENSION; c++)
    if (s->u[c])
      gts_object_destroy (GTS_OBJECT (s->u[c]));

  (* GTS_OBJECT_CLASS (gfs_source_class ())->parent_class->destroy) (o);
}

static void gfs_source_coriolis_read (GtsObject ** o, GtsFile * fp)
{
  FttComponent c;
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

  (* GTS_OBJECT_CLASS (gfs_source_coriolis_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  for (c = 0; c < FTT_DIMENSION; c++) {
    GfsVariable * v = GFS_SOURCE_VELOCITY (*o)->v[c];

    if (v->sources) {
      GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;
      
      while (i) {
	if (i->data != *o && GFS_IS_SOURCE_CORIOLIS (i->data)) {
	  gts_file_error (fp, "variable '%s' cannot have multiple Coriolis source terms", v->name);
	  return;
	}
	i = i->next;
      }
    }
  }

  GfsSourceCoriolis * s = GFS_SOURCE_CORIOLIS (*o);

  s->omegaz = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_read (s->omegaz, gfs_object_simulation (s), fp);

  if (fp->type != '\n') {
    s->drag = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_read (s->drag, gfs_object_simulation (s), fp);
  }

  if (s->beta < 1.)
    for (c = 0; c <  FTT_DIMENSION; c++)
      s->u[c] = gfs_temporary_variable (domain);
  else {
    GFS_SOURCE_GENERIC (s)->centered_value = NULL;
    GFS_SOURCE_GENERIC (s)->mac_value = NULL;
  }

#if !FTT_2D
  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_DOUBLE, "x", TRUE, &s->d[0]},
      {GTS_DOUBLE, "y", TRUE, &s->d[1]},
      {GTS_DOUBLE, "z", TRUE, &s->d[2]},
      {GTS_NONE}
    };
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
    gdouble norm = sqrt (s->d[0]*s->d[0] + s->d[1]*s->d[1] + s->d[2]*s->d[2]);
    if (norm == 0.)
      gts_file_error (fp, "vector must be non-zero");
    else {
      s->d[0] /= norm;
      s->d[1] /= norm;
      s->d[2] /= norm;
    }      
  }
#endif /* 3D */
}

static void gfs_source_coriolis_write (GtsObject * o, FILE * fp)
{
  GfsSourceCoriolis * s = GFS_SOURCE_CORIOLIS (o);

  (* GTS_OBJECT_CLASS (gfs_source_coriolis_class ())->parent_class->write) (o, fp);
  gfs_function_write (s->omegaz, fp);
  if (s->drag)
    gfs_function_write (s->drag, fp);
  if (s->d[2] != 1.)
    fprintf (fp," { x = %g y = %g z = %g }", s->d[0], s->d[1], s->d[2]);
}

static void coriolis_drag_metric_coefficients (GfsSourceCoriolis * sc,
					       FttCell * cell,
					       GfsVariable ** u,
					       FttComponent c1,
					       gdouble f[2])
{
  f[0] = sc->drag ? gfs_function_value (sc->drag, cell) : 0.;
  f[1] = sc->omegaz ? gfs_function_value (sc->omegaz, cell) : 0.;

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (sc));
  if (domain->advection_metric) {
    /* fixme: 2D only */
    FttComponent c2 = (c1 + 1) % 2;
    gdouble m[2];
    (* domain->advection_metric) (domain, cell, c1, m);
    f[0] += GFS_VALUE (cell, u[c2])*m[0];
    f[1] += (c1 == FTT_X ? 1. : -1.)*GFS_VALUE (cell, u[c2])*m[1];
  }
}

static gdouble gfs_source_coriolis_mac_value (GfsSourceGeneric * s,
					      FttCell * cell,
					      GfsVariable * v)
{
  GfsSourceVelocity * sv = GFS_SOURCE_VELOCITY (s);
  GfsSourceCoriolis * sc = GFS_SOURCE_CORIOLIS (s);

#if FTT_2D
  gdouble f[2];
  coriolis_drag_metric_coefficients (sc, cell, sv->v, v->component, f);
  switch (v->component) {
  case FTT_X: return   f[1]*GFS_VALUE (cell, sv->v[1]) - f[0]*GFS_VALUE (cell, sv->v[0]);
  case FTT_Y: return - f[0]*GFS_VALUE (cell, sv->v[1]) - f[1]*GFS_VALUE (cell, sv->v[0]);
  default: g_assert_not_reached ();
  }
#else  /* 3D */
  gdouble f = gfs_function_value (sc->omegaz, cell);
  gdouble e = sc->drag ? gfs_function_value (sc->drag, cell) : 0.;

  switch (v->component) {
  case FTT_X: return - e*GFS_VALUE (cell, sv->v[0]) +
      f*(GFS_VALUE (cell, sv->v[1])*sc->d[2] - GFS_VALUE (cell, sv->v[2])*sc->d[1]);
  case FTT_Y: return - e*GFS_VALUE (cell, sv->v[1]) +
      f*(GFS_VALUE (cell, sv->v[2])*sc->d[0] - GFS_VALUE (cell, sv->v[0])*sc->d[2]);
  case FTT_Z: return - e*GFS_VALUE (cell, sv->v[2]) +
      f*(GFS_VALUE (cell, sv->v[0])*sc->d[1] - GFS_VALUE (cell, sv->v[1])*sc->d[0]);
  default: g_assert_not_reached ();
  }
#endif /* 3D */
  return 0.;
}

static void save_coriolis (FttCell * cell, GfsSourceCoriolis * s)
{
  GfsSourceVelocity * sv = GFS_SOURCE_VELOCITY (s);
#if FTT_2D
  gdouble f[2];
  coriolis_drag_metric_coefficients (s, cell, sv->v, FTT_X, f);
  GFS_VALUE (cell, s->u[0]) =   
    (1. - s->beta)*(f[1]*GFS_VALUE (cell, sv->v[1]) - f[0]*GFS_VALUE (cell, sv->v[0]));
  coriolis_drag_metric_coefficients (s, cell, sv->v, FTT_Y, f);
  GFS_VALUE (cell, s->u[1]) = 
    (1. - s->beta)*(- f[0]*GFS_VALUE (cell, sv->v[1]) - f[1]*GFS_VALUE (cell, sv->v[0]));
#else  /* 3D */
  gdouble f = gfs_function_value (s->omegaz, cell)*(1. - s->beta);
  gdouble e = s->drag ? gfs_function_value (s->drag, cell)*(1. - s->beta) : 0.;

  GFS_VALUE (cell, s->u[0]) = - e*GFS_VALUE (cell, sv->v[0]) +
    f*(GFS_VALUE (cell, sv->v[1])*s->d[2] - GFS_VALUE (cell, sv->v[2])*s->d[1]);
  GFS_VALUE (cell, s->u[1]) = - e*GFS_VALUE (cell, sv->v[1]) +
    f*(GFS_VALUE (cell, sv->v[2])*s->d[0] - GFS_VALUE (cell, sv->v[0])*s->d[2]);
  GFS_VALUE (cell, s->u[2]) = - e*GFS_VALUE (cell, sv->v[2]) +
    f*(GFS_VALUE (cell, sv->v[0])*s->d[1] - GFS_VALUE (cell, sv->v[1])*s->d[0]);
#endif /* 3D */
}

static gboolean gfs_source_coriolis_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class)->event) 
      (event, sim)) {
    if (GFS_SOURCE_CORIOLIS (event)->beta < 1.) {
      gfs_catch_floating_point_exceptions ();
      gfs_domain_traverse_layers (GFS_DOMAIN (sim), (FttCellTraverseFunc) save_coriolis, event);
      if (gfs_restore_floating_point_exceptions ()) {
	GfsSourceCoriolis * c = GFS_SOURCE_CORIOLIS (event);
	gchar * s = g_strconcat ("\n", gfs_function_description (c->omegaz, FALSE), NULL);
	if (c->drag)
	  s = g_strconcat (s, "\n", gfs_function_description (c->drag, FALSE), NULL);
	/* fixme: memory leaks */
	g_message ("floating-point exception in user-defined function(s):%s", s);
	exit (1);
      }
    }
    return TRUE;
  }
  return FALSE;
}

static gdouble gfs_source_coriolis_centered_value (GfsSourceGeneric * s,
						   FttCell * cell,
						   GfsVariable * v)
{
  return GFS_VALUE (cell, GFS_SOURCE_CORIOLIS (s)->u[v->component]);
}

static void gfs_source_coriolis_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = source_coriolis_destroy;
  GTS_OBJECT_CLASS (klass)->read = gfs_source_coriolis_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_source_coriolis_write;

  GFS_EVENT_CLASS (klass)->event = gfs_source_coriolis_event;
}

static void gfs_source_coriolis_init (GfsSourceGeneric * s)
{
  s->mac_value = gfs_source_coriolis_mac_value;
  GFS_SOURCE_CORIOLIS (s)->beta = 0.5; /* Crank-Nicholson */
  s->centered_value = gfs_source_coriolis_centered_value;
  GfsSourceCoriolis * sc = GFS_SOURCE_CORIOLIS (s);
  sc->d[0] = sc->d[1] = 0.; sc->d[2] = 1.;
}

GfsSourceGenericClass * gfs_source_coriolis_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_coriolis_info = {
      "GfsSourceCoriolis",
      sizeof (GfsSourceCoriolis),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_coriolis_class_init,
      (GtsObjectInitFunc) gfs_source_coriolis_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_velocity_class ()),
				  &gfs_source_coriolis_info);
  }

  return klass;
}

/**
 * gfs_has_source_coriolis:
 * @domain: a #GfsDomain.
 *
 * Returns: the #GfsSourceCoriolis associated with @domain or %NULL.
 */
GfsSourceCoriolis * gfs_has_source_coriolis (GfsDomain * domain)
{
  GfsVariable * v;

  g_return_val_if_fail (domain != NULL, NULL);

  v = gfs_variable_from_name (domain->variables, "U");
  g_return_val_if_fail (v != NULL, NULL);

  if (v->sources) {
    GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;

    while (i) {
      if (GFS_IS_SOURCE_CORIOLIS (i->data))
	return i->data;
      i = i->next;
    }
  }
  return NULL;
}

static void implicit_coriolis_2D (FttCell * cell, GfsSourceCoriolis * s)
{
  GfsSourceVelocity * sv = GFS_SOURCE_VELOCITY (s);
  gdouble m[2][2];
  gdouble f[2];
  GfsSimulation * sim = gfs_object_simulation (s);
  gdouble dt = sim->advection_params.dt*s->beta;

  coriolis_drag_metric_coefficients (s, cell, sv->v, FTT_X, f);
  m[0][0] = 1. + f[0]*dt;
  m[0][1] = - f[1]*dt;
  coriolis_drag_metric_coefficients (s, cell, sv->v, FTT_Y, f);
  m[1][0] =   f[1]*dt;
  m[1][1] = 1. + f[0]*dt;

  gdouble det = m[0][0]*m[1][1] - m[0][1]*m[1][0];
  gdouble u = GFS_VALUE (cell, sv->v[0]);
  gdouble v = GFS_VALUE (cell, sv->v[1]);

  GFS_VALUE (cell, sv->v[0]) = (  m[1][1]*u - m[0][1]*v)/det;
  GFS_VALUE (cell, sv->v[1]) = (- m[1][0]*u + m[0][0]*v)/det;
}

static void implicit_coriolis_3D (FttCell * cell, GfsSourceCoriolis * s)
{
  GfsSourceVelocity * sv = GFS_SOURCE_VELOCITY (s);
  GfsSimulation * sim = gfs_object_simulation (s);
  gdouble c = sim->advection_params.dt*gfs_function_value (s->omegaz, cell)*s->beta;
  gdouble e = s->drag ? sim->advection_params.dt*gfs_function_value (s->drag, cell)*s->beta : 0.;
  GtsMatrix * m = gts_matrix_new (    1. + e, -c*s->d[2],   c*s->d[1],  0.,
				   c*s->d[2],     1. + e,  -c*s->d[0],  0.,
				  -c*s->d[1],  c*s->d[0],      1. + e,  0.,
				          0.,         0.,          0.,  0.);
  GtsMatrix * im = gts_matrix3_inverse (m);
  g_assert (im);
  gint i, j;
  gdouble vel[3];

  for (i = 0; i < 3; i++)
    vel[i] = GFS_VALUE (cell, sv->v[i]);
  
  for (i = 0; i < 3; i++) {
    GFS_VALUE (cell, sv->v[i]) = 0.;
    for (j = 0; j < 3; j++)
      GFS_VALUE (cell, sv->v[i]) += im[i][j]*vel[j];
  }

  gts_matrix_destroy (m);
  gts_matrix_destroy (im);
}

/**
 * gfs_source_coriolis_implicit:
 * @domain: a #GfsDomain.
 * @dt: the timestep.
 *
 * Applies the implicit part of the Coriolis source term of @domain.
 */
void gfs_source_coriolis_implicit (GfsDomain * domain,
				   gdouble dt)
{
  GfsSourceCoriolis * s;

  g_return_if_fail (domain != NULL);

  if ((s = gfs_has_source_coriolis (domain))) {
    GfsSimulation * sim = GFS_SIMULATION (domain);
    gdouble olddt = sim->advection_params.dt;
    sim->advection_params.dt = dt;
    gfs_catch_floating_point_exceptions ();
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) (s->d[2] != 1. ? 
						     implicit_coriolis_3D : implicit_coriolis_2D),
			      s);
    if (gfs_restore_floating_point_exceptions ()) {
      gchar * c = g_strconcat ("\n", gfs_function_description (s->omegaz, FALSE), NULL);
      if (s->drag)
	c = g_strconcat (c, "\n", gfs_function_description (s->drag, FALSE), NULL);
      /* fixme: memory leaks */
      g_message ("floating-point exception in user-defined function(s):%s", c);
      exit (1);
    }
    sim->advection_params.dt = olddt;
  }
}

/** \endobject{GfsSourceCoriolis} */
