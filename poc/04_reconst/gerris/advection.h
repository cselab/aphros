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

#ifndef __ADVECTION_H__
#define __ADVECTION_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "domain.h"
#include "poisson.h"

#define GFS_SMALL 0.5

typedef enum {
  GFS_GODUNOV,
  GFS_NONE
} GfsAdvectionScheme;

typedef enum {
  GFS_CENTERED_UPWINDING,
  GFS_FACE_UPWINDING,
  GFS_NO_UPWINDING,
} GfsUpwinding;

typedef struct _GfsAdvectionParams GfsAdvectionParams;
typedef
void      (* GfsFaceAdvectionFluxFunc)       (const FttCellFace * face,
					      const GfsAdvectionParams * par);
typedef void (* GfsMergedTraverseFunc)       (GSList * merged,
					      gpointer data);

struct _GfsAdvectionParams {
  gdouble cfl, dt;
  GfsVariable * v, * fv, ** u, ** g;
  GfsCenterGradient gradient;
  gboolean use_centered_velocity;
  GfsUpwinding upwinding;
  GfsFaceAdvectionFluxFunc flux;
  GfsAdvectionScheme scheme;
  gboolean average, gc;
  GfsMergedTraverseFunc update;
  guint moving_order;
  GfsFunction * sink[FTT_DIMENSION];
  gboolean linear;
  void (* diffusion_solve) (GfsDomain * domain,
			    GfsMultilevelParams * par,
			    GfsVariable * v,
			    GfsVariable * rhs, 
			    GfsVariable * rhoc,
			    GfsVariable * axi);
};

void         gfs_advection_params_init        (GfsAdvectionParams * par);
void         gfs_advection_params_write       (GfsAdvectionParams * par, 
					       FILE * fp);
void         gfs_advection_params_read        (GfsAdvectionParams * par, 
					       GtsFile * fp);
void         gfs_cell_advected_face_values    (FttCell * cell,
					       const GfsAdvectionParams * par);
void         gfs_cell_non_advected_face_values (FttCell * cell,
						const GfsAdvectionParams * par);
gdouble      gfs_face_upwinded_value          (const FttCellFace * face,
					       GfsUpwinding upwinding,
					       GfsVariable ** u);
void         gfs_face_advection_flux          (const FttCellFace * face,
					       const GfsAdvectionParams * par);
void         gfs_face_velocity_advection_flux (const FttCellFace * face,
					       const GfsAdvectionParams * par);
void         gfs_face_velocity_convective_flux (const FttCellFace * face,
						const GfsAdvectionParams * par);
void         gfs_face_advected_normal_velocity     (const FttCellFace * face,
						    const GfsAdvectionParams * par);
void         gfs_face_interpolated_normal_velocity (const FttCellFace * face,
						    GfsVariable ** v);
void         gfs_face_reset_normal_velocity        (const FttCellFace * face);
gboolean     gfs_cell_is_small                     (const FttCell * cell);
void         gfs_set_merged                        (GfsDomain * domain);
void         gfs_domain_traverse_merged            (GfsDomain * domain,
						    GfsMergedTraverseFunc func,
						    gpointer data);
void         gfs_advection_update                  (GSList * merged, 
					            const GfsAdvectionParams * par);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __ADVECTION_H__ */
