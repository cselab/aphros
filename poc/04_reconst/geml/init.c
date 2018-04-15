/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2012 National Institute of Water and Atmospheric Research
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
/** \file
 * \brief Initialisation.
 */

/** \mainpage Programming API Reference
 *
 * - <a href=modules.html>Class hierarchy</a>
 * - <a href=globals.html>Function index</a>
 * - <a href=classes.html>Data structure index</a>
 */

#include "config.h"

#ifdef HAVE_FENV_H
# define _GNU_SOURCE
# include <fenv.h>
# ifdef FE_NOMASK_ENV
#  ifdef FE_DIVBYZERO
#    ifdef FE_INVALID
#     define EXCEPTIONS (FE_DIVBYZERO|FE_INVALID)
#    else
#     define EXCEPTIONS (FE_DIVBYZERO)
#   endif
#  endif /* !FE_DIVBYZERO */
# endif /* FE_NO_MASK_ENV */
#endif /* HAVE_FENV_H */

#include <stdlib.h>
#include <locale.h>

#include "boundary.h"
#include "mpi_boundary.h"
#include "init.h"
#include "refine.h"
#include "output.h"
#include "adaptive.h"
#include "source.h"
#include "tension.h"
#include "ocean.h"
#include "wave.h"
#include "levelset.h"
#include "vof.h"
#include "solid.h"
#include "moving.h"
#include "river.h"
#include "balance.h"
#include "map.h"
#include "metric.h"
#include "particle.h"
#include "cartesian.h"

#ifdef HAVE_MPI
# include <mpi.h>
#endif /* HAVE_MPI */

static void gfs_log (const gchar * log_domain,
		     GLogLevelFlags log_level,
		     const gchar * message)
{
  int type = 0;
  gchar * pe;
  const gchar stype[][10] = {
    "ERROR", "CRITICAL", "WARNING", "MESSAGE", "INFO", "DEBUG"
  };

#ifdef HAVE_MPI
  int rank = -1;
  MPI_Comm_size (MPI_COMM_WORLD, &rank);
  if (rank > 1)
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  else
    rank = -1;
  if (rank >= 0) {
    char name[MPI_MAX_PROCESSOR_NAME];
    int length;
    MPI_Get_processor_name (name, &length);
    pe = g_strdup_printf ("PE %d (%s): ", rank, name);
  }
  else
#endif /* HAVE_MPI */
    pe = g_strdup ("");

  switch (log_level & G_LOG_LEVEL_MASK) {
  case G_LOG_LEVEL_ERROR:    type = 0; break;
  case G_LOG_LEVEL_CRITICAL: type = 1; break;
  case G_LOG_LEVEL_WARNING:  type = 2; break;
  case G_LOG_LEVEL_MESSAGE:  type = 3; break;
  case G_LOG_LEVEL_INFO:     type = 4; break;
  case G_LOG_LEVEL_DEBUG:    type = 5; break;
  default:
    g_assert_not_reached ();
  }
  fprintf (stderr, "\n%s-%s **: %s%s\n\n", 
	   log_domain, stype[type], pe, message);
  g_free (pe);
}

/**
 * gfs_classes:
 *
 * Returns: a pointer to a NULL-terminated array of all the classes
 * usable in Gerris parameter files.
 */
GtsObjectClass ** gfs_classes (void)
{
  static GtsObjectClass ** classes = NULL;
  if (classes == NULL) { gpointer klass[] = {

  gfs_global_class (),
  gfs_domain_class (),
    gfs_simulation_class (),
      gfs_ocean_class (),
      gfs_advection_class (),
      gfs_poisson_class (),
      gfs_simulation_moving_class (),
      gfs_axi_class (),
        gfs_advection_axi_class (),
      gfs_wave_class (),
      gfs_river_class (),
    gfs_domain_projection_class (),

  gfs_surface_bc_class (),

  gfs_box_class (),

  gfs_gedge_class (),

  gfs_bc_dirichlet_class (),
  gfs_bc_subcritical_class (),
  gfs_bc_neumann_class (),
    gfs_bc_angle_class (),
  gfs_bc_navier_class (),
  gfs_bc_flather_class (),

  gfs_boundary_class (),
    gfs_boundary_inflow_constant_class (),
    gfs_boundary_outflow_class (),
    gfs_boundary_gradient_class (),
    gfs_boundary_periodic_class (),
      gfs_boundary_mpi_class (),

  gfs_refine_class (),
    gfs_refine_solid_class (),
    gfs_refine_surface_class (),
      gfs_refine_distance_class (),
      gfs_refine_height_class (),
    gfs_layers_class (),

  gfs_event_class (),
    gfs_variable_class (),
      gfs_variable_boolean_class (),
      gfs_variable_tracer_class (),
        gfs_variable_vof_concentration_class (),
        gfs_variable_tracer_vof_class (),
          gfs_variable_tracer_vof_height_class (),
      gfs_variable_residual_class (),
      gfs_variable_filtered_class (),
      gfs_variable_diagonal_class (),
      gfs_variable_function_class (),
#if FTT_2D
        gfs_variable_stream_function_class (),
#endif /* FTT_2D */
        gfs_variable_average_class (),
        gfs_variable_poisson_class (),
        gfs_variable_laplacian_class (),
      gfs_hydrostatic_pressure_class (),
      gfs_variable_age_class (),
      gfs_variable_curvature_class (),
        gfs_variable_position_class (),
      gfs_variable_distance_class (),

    gfs_constant_class (),
      gfs_discharge_elevation_class (),
      gfs_spatial_sum_class (),

    gfs_solid_class (),
      gfs_solid_moving_class(),

    gfs_init_class (),
    gfs_init_mask_class (),
    gfs_init_flow_constant_class (),
    gfs_init_fraction_class (),
    gfs_init_vorticity_class (),
    gfs_init_wave_class (),

    gfs_generic_metric_class (),
      gfs_metric_stretch_class (),
      gfs_variable_metric_class (),
        gfs_metric_lon_lat_class (),
        gfs_stored_metric_class (),
          gfs_metric_class (),
          gfs_metric_cubed_class (),
          gfs_metric_cubed1_class (),
          gfs_metric_variable_class (),
            gfs_metric_laplace_class (),

    gfs_adapt_class (),
      gfs_adapt_vorticity_class (),
      gfs_adapt_streamline_curvature_class (),
      gfs_adapt_function_class (),
      gfs_adapt_thickness_class (),
      gfs_adapt_gradient_class (),
        gfs_adapt_error_class (),

    gfs_event_sum_class (),
      gfs_event_sum_direction_class (),
    gfs_event_harmonic_class (),
    gfs_event_stop_class (),
    gfs_event_script_class (),
    gfs_event_balance_class (),
    gfs_source_generic_class (),
      gfs_source_scalar_class (),
        gfs_source_class (),
          gfs_source_control_class (),
            gfs_source_control_field_class (),
          gfs_source_flux_class (),
          gfs_source_pipe_class (),
        gfs_source_diffusion_class (),
          gfs_source_diffusion_explicit_class (),
      gfs_source_velocity_class (),
        gfs_source_viscosity_class (),
          gfs_source_viscosity_explicit_class (),
        gfs_source_friction_class (),
        gfs_source_coriolis_class (),
          gfs_source_tension_class (),
          gfs_source_tension_css_class (),
#if !FTT_2D
        gfs_source_hydrostatic_class (),
#endif /* 3D */
    gfs_remove_droplets_class (),
    gfs_remove_ponds_class (),
    gfs_event_filter_class (),
    gfs_event_list_class (),

    gfs_diffusion_class (),
   
    gfs_output_class (),
      gfs_output_time_class (),
      gfs_output_progress_class (),
      gfs_output_projection_stats_class (),
      gfs_output_diffusion_stats_class (),
      gfs_output_solid_stats_class (),
      gfs_output_adapt_stats_class (),
      gfs_output_timing_class (),
      gfs_output_balance_class (),
      gfs_output_solid_force_class (),
      gfs_output_location_class (),
        gfs_output_particle_class (),
      gfs_output_simulation_class (),
      gfs_output_boundaries_class (),
      gfs_output_object_class (),

      gfs_output_scalar_class (),
        gfs_output_scalar_norm_class (),
        gfs_output_scalar_stats_class (),
        gfs_output_scalar_sum_class (),
        gfs_output_scalar_maxima_class (),
        gfs_output_scalar_histogram_class (),
        gfs_output_droplet_sums_class (),
        gfs_output_error_norm_class (),
          gfs_output_correlation_class (),
	gfs_output_squares_class (),
	gfs_output_streamline_class (),
        gfs_output_ppm_class (),  
        gfs_output_grd_class (),  

  gfs_map_class (),
    gfs_map_function_class (),
    gfs_map_transform_class (),

  gfs_particle_class (),

  gfs_cartesian_grid_class (),

  gfs_derived_variable_class (),

  gfs_generic_surface_class (),
    gfs_surface_class (),

  gfs_function_class (),
    gfs_function_constant_class (),
    gfs_function_spatial_class (),
      gfs_function_map_class (),

  NULL};

    guint n = 0;
    gpointer * c = klass;
    while (*(c++)) n++;
    classes = g_malloc ((n + 1)*sizeof (gpointer));
    memcpy (classes, klass, (n + 1)*sizeof (gpointer));
  }
  return classes;
}

static gboolean disabled_fpe = FALSE;

typedef void (* AtExitFunc) (void);

/**
 * gfs_catch_floating_point_exceptions:
 *
 * Catch the default floating-point exceptions set in the Gerris
 * library.
 */
void gfs_catch_floating_point_exceptions (void)
{
#ifdef EXCEPTIONS
  fedisableexcept (EXCEPTIONS);
  feclearexcept (EXCEPTIONS);
#endif /* EXCEPTIONS */
}

/**
 * gfs_restore_floating_point_exceptions:
 *
 * Restores the default floating-point exceptions set in the Gerris
 * library.
 *
 * Returns: 0 if no exceptions where raised between this call and the
 * call to gfs_catch_floating_point_exceptions(), non-zero otherwise.
 */
int gfs_restore_floating_point_exceptions (void)
{
#ifdef EXCEPTIONS
  int ret = fetestexcept (EXCEPTIONS);
  feclearexcept (EXCEPTIONS);
  if (!disabled_fpe)
    feenableexcept (EXCEPTIONS);
  return ret;
#else /* !EXCEPTIONS */
  return 0;
#endif /* !EXCEPTIONS */
}

/**
 * gfs_disable_floating_point_exceptions:
 *
 * Disables floating-point exceptions (they are enabled by default
 * when using the Gerris library).
 */
void gfs_disable_floating_point_exceptions (void)
{
#ifdef EXCEPTIONS
  disabled_fpe = TRUE;
  fedisableexcept (EXCEPTIONS);
#endif /* !EXCEPTIONS */
}

/**
 * gfs_enable_floating_point_exceptions:
 *
 * Enables floating-point exceptions (they are enabled by default
 * when using the Gerris library).
 */
void gfs_enable_floating_point_exceptions (void)
{
#ifdef EXCEPTIONS
  disabled_fpe = FALSE;
  feenableexcept (EXCEPTIONS);
#endif /* !EXCEPTIONS */
}

/**
 * gfs_init:
 * @argc: a pointer on the number of command line arguments passed to
 * the program.
 * @argv: a pointer on the command line arguments passed to the
 * program.
 *
 * Initializes the Gerris library. This function must be called before
 * any other function of the library.
 */
void gfs_init (int * argc, char *** argv)
{
  static gboolean initialized = FALSE;

  if (initialized)
    return;

  if (!setlocale (LC_ALL, "POSIX"))
    g_warning ("cannot set locale to POSIX");

#ifdef HAVE_MPI
  MPI_Initialized (&initialized);
  if (!initialized) {
    if (!argc || !argv) {
      int argc1 = 1;
      char ** argv1;

      argv1 = g_malloc (sizeof (char *));
      argv1[0] = g_strdup ("gfs_init");
      MPI_Init (&argc1, &argv1);
      g_free (argv1[0]); g_free (argv1);
    }
    else
      MPI_Init (argc, argv);
    MPI_Errhandler_set (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
  }
  atexit ((AtExitFunc) MPI_Finalize);
#endif /* HAVE_MPI */
  initialized = TRUE;

#ifdef EXCEPTIONS
  feenableexcept (EXCEPTIONS);
#endif /* EXCEPTIONS */

  g_log_set_handler (G_LOG_DOMAIN,
		     G_LOG_LEVEL_ERROR |
		     G_LOG_LEVEL_CRITICAL |
		     G_LOG_LEVEL_WARNING |
		     G_LOG_LEVEL_MESSAGE |
		     G_LOG_LEVEL_INFO |
		     G_LOG_LEVEL_DEBUG |
		     G_LOG_FLAG_FATAL |
		     G_LOG_FLAG_RECURSION,
		     (GLogFunc) gfs_log, NULL);

  /* Instantiates classes before reading any domain or simulation file */
  gfs_classes ();
}
