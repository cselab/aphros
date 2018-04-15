#ifndef __SOLID_H__
#define __SOLID_H__

#include <gts.h>

#include "domain.h"
#include "event.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

enum {
  GFS_STATUS_UNDEFINED = 0,
  GFS_STATUS_SOLID     = 1,
  GFS_STATUS_FLUID     = 2
};

void         gfs_cell_fluid                              (FttCell * cell);
gboolean     gfs_solid_is_thin                           (FttCell * cell, 
							  GfsGenericSurface * s);
gboolean     gfs_set_2D_solid_fractions_from_surface     (FttCell * cell,
							  GfsGenericSurface * s);
guint        gfs_init_solid_fractions_leaves             (GfsDomain * domain,
							  GSList * i,
							  GfsVariable * status);
void         gfs_init_solid_fractions_from_children      (GfsDomain * domain,
							  gboolean destroy_solid,
							  FttCellCleanupFunc cleanup,
							  gpointer data,
							  GfsVariable * status);
guint        gfs_domain_init_solid_fractions             (GfsDomain * domain,
							  GSList * i,
							  gboolean destroy_solid,
							  FttCellCleanupFunc cleanup,
							  gpointer data,
							  GfsVariable * status);
void         gfs_cell_init_solid_fractions_from_children (FttCell * cell);
gboolean     gfs_cell_check_solid_fractions              (FttCell * root);
void         gfs_domain_init_fraction                    (GfsDomain * domain,
							  GfsGenericSurface * s,
							  GfsVariable * c);
void         gfs_cell_cm                                 (const FttCell * cell, 
							  FttVector * cm);
void         gfs_solid_normal                            (const FttCell * cell,
							  FttVector * n);
void         gfs_face_ca                                 (const FttCellFace * face, 
							  FttVector * ca);
void         gfs_solid_coarse_fine                       (FttCell * parent,
							  GfsDomain * domain);

/* GfsSolid: Header */

typedef struct _GfsSolid         GfsSolid;

struct _GfsSolid {
  /*< private >*/
  GfsEvent parent;

  /*< public >*/
  GfsGenericSurface * s;
};

#define GFS_SOLID(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSolid,\
					         gfs_solid_class ())
#define GFS_IS_SOLID(obj)         (gts_object_is_from_class (obj,\
						 gfs_solid_class ()))

GfsEventClass * gfs_solid_class  (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __SOLID_H__ */
