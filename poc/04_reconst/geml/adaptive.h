#ifndef __ADAPTIVE_H__
#define __ADAPTIVE_H__

#include "simulation.h"
#include "event.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void          gfs_cell_coarse_init          (FttCell * cell,
					     GfsDomain * domain);
void          gfs_adapt_stats_init          (GfsAdaptStats * s);
void          gfs_adapt_stats_update        (GfsAdaptStats * s);
gboolean      gfs_simulation_adapt          (GfsSimulation * simulation);
void          gfs_domain_reshape            (GfsDomain * domain,
					     guint depth);

/* GfsAdapt: Header */

typedef struct _GfsAdapt         GfsAdapt;

struct _GfsAdapt {
  /*< private >*/
  GfsEvent parent;
  gboolean active;

  /*< public >*/
  GfsFunction * minlevel, * maxlevel;
  guint mincells, maxcells;
  gdouble cmax, weight, cfactor;
  GfsVariable * c;
  GtsKeyFunc cost;
};

#define GFS_ADAPT(obj)            GTS_OBJECT_CAST (obj,\
					         GfsAdapt,\
					         gfs_adapt_class ())
#define GFS_IS_ADAPT(obj)         (gts_object_is_from_class (obj,\
						 gfs_adapt_class ()))

GfsEventClass * gfs_adapt_class  (void);

/* GfsAdaptVorticity: Header */

typedef struct _GfsAdaptVorticity         GfsAdaptVorticity;

struct _GfsAdaptVorticity {
  /*< private >*/
  GfsAdapt parent;
  GfsVariable ** u;
  gdouble maxa;

  /*< public >*/
};

#define GFS_ADAPT_VORTICITY(obj)            GTS_OBJECT_CAST (obj,\
					         GfsAdaptVorticity,\
					         gfs_adapt_vorticity_class ())
#define GFS_IS_ADAPT_VORTICITY(obj)         (gts_object_is_from_class (obj,\
						 gfs_adapt_vorticity_class ()))

GfsEventClass * gfs_adapt_vorticity_class  (void);
 
/* GfsAdaptStreamlineCurvature: Header */

#define GFS_IS_ADAPT_STREAMLINE_CURVATURE(obj)         (gts_object_is_from_class (obj,\
						 gfs_adapt_streamline_curvature_class ()))

GfsEventClass * gfs_adapt_streamline_curvature_class  (void);
 
/* GfsAdaptFunction: Header */

typedef struct _GfsAdaptFunction         GfsAdaptFunction;

struct _GfsAdaptFunction {
  /*< private >*/
  GfsAdapt parent;

  /*< public >*/
  GfsFunction * f;
};

#define GFS_ADAPT_FUNCTION(obj)            GTS_OBJECT_CAST (obj,\
					         GfsAdaptFunction,\
					         gfs_adapt_function_class ())
#define GFS_IS_ADAPT_FUNCTION(obj)         (gts_object_is_from_class (obj,\
						 gfs_adapt_function_class ()))

GfsEventClass * gfs_adapt_function_class  (void);

/* GfsAdaptGradient: Header */

typedef struct _GfsAdaptGradient         GfsAdaptGradient;

struct _GfsAdaptGradient {
  /*< private >*/
  GfsAdaptFunction parent;
  gdouble dimension;

  /*< public >*/
  GfsVariable * v;
};

#define GFS_ADAPT_GRADIENT(obj)            GTS_OBJECT_CAST (obj,\
					         GfsAdaptGradient,\
					         gfs_adapt_gradient_class ())
#define GFS_IS_ADAPT_GRADIENT(obj)         (gts_object_is_from_class (obj,\
						 gfs_adapt_gradient_class ()))

GfsEventClass * gfs_adapt_gradient_class  (void);

/* GfsAdaptError: Header */

typedef struct _GfsAdaptError         GfsAdaptError;

struct _GfsAdaptError {
  /*< private >*/
  GfsAdaptGradient parent;
  GfsVariable * dv[FTT_DIMENSION];
  FttComponent c;

  /*< public >*/
  GfsVariable * v;
};

#define GFS_ADAPT_ERROR(obj)            GTS_OBJECT_CAST (obj,\
					         GfsAdaptError,\
					         gfs_adapt_error_class ())
#define GFS_IS_ADAPT_ERROR(obj)         (gts_object_is_from_class (obj,\
						 gfs_adapt_error_class ()))

GfsEventClass * gfs_adapt_error_class  (void);

/* GfsAdaptThickness: Header */

typedef struct _GfsAdaptThickness         GfsAdaptThickness;

struct _GfsAdaptThickness {
  /*< private >*/
  GfsAdapt parent;

  /*< public >*/
  GfsVariable * v, * c;
};

#define GFS_ADAPT_THICKNESS(obj)            GTS_OBJECT_CAST (obj,\
					         GfsAdaptThickness,\
					         gfs_adapt_thickness_class ())
#define GFS_IS_ADAPT_THICKNESS(obj)         (gts_object_is_from_class (obj,\
						 gfs_adapt_thickness_class ()))

GfsEventClass * gfs_adapt_thickness_class  (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __ADAPTIVE_H__ */
