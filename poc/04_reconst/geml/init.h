#ifndef __INIT_H__
#define __INIT_H__

#include <gts.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

GtsObjectClass ** gfs_classes             (void);
void              gfs_init                (int * argc, 
					   char *** argv);
void gfs_catch_floating_point_exceptions   (void);
int  gfs_restore_floating_point_exceptions (void);
void gfs_disable_floating_point_exceptions (void);
void gfs_enable_floating_point_exceptions  (void);

#define gfs_restore_fpe_for_function(f) \
       { \
         if (gfs_restore_floating_point_exceptions ()) { \
           g_message ("floating-point exception in user-defined function:\n%s", \
	              gfs_function_description (f, FALSE)); \
           exit (1); \
         } \
       }

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __INIT_H__ */
