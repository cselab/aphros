#include "grid/octree.h"
#include "utils.h"
#include "distance.h"
#include "fractions.h"
#include "view.h"
#include <inside.h>

static const char *me = "distance";
static void
usg(void) {
  fprintf(stderr, "%s [off|ply] < off\n", me);
  exit(1);
}

int main(int argc, char **argv)
{
  coord max;
  coord min;
  coord * p;
  double maxl;
  double *ver;
  int i;
  int nt;
  int nv;
  int *tri;

  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
        usg();
        break;
    default:
        fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
        exit(2);
    }
  if (argv[0] == NULL) {
      fprintf(stderr, "%s: needs a mesh file\n", me);
      exit(2);
  }
  if (inside_mesh_read(argv[0], &nt, &tri, &nv, &ver) != 0) {
      fprintf(stderr, "%s: fail to read '%s'\n", me, argv[0]);
      exit(2);
  }
  if ((p = malloc( (3 * nt + 1) * sizeof *p)) == NULL) {
      fprintf(stderr, "%s: malloc failed\n", me);
      exit(2);
  }
  for (i = 0; i < 3 * nt; i++) {
      p[i].x = ver[3 * tri[i]];
      p[i].y = ver[3 * tri[i] + 1];
      p[i].z = ver[3 * tri[i] + 2];
  }
  p[3 * nt].x = nodata;
  free(tri);
  free(ver);
  maxl = -HUGE;
  bounding_box (p, &min, &max);
  foreach_dimension()
    if (max.x - min.x > maxl)
      maxl = max.x - min.x;
  init_grid (4);
  size (1.2*maxl);
  origin ((max.x + min.x)/2. - L0/2,
          (max.y + min.y)/2. - L0/2,
          (max.z + min.z)/2. - L0/2);

  scalar d[];
  distance (d, p);
  while (adapt_wavelet ({d}, (double[]){5e-4*L0}, 10).nf);

  view (fov = 15.65, quat = {-0.52,0.31,0.38,-0.7},
        tx = -0.045, ty = 0.015, width = 640, height = 480, bg = {1,1,1});
  isosurface ("d", 0, color = "level", min = 5, max = 10);
  save ("level.png");

  vertex scalar phi[];
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
             d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  boundary ({phi});
  scalar f[];
  face vector s[];
  fractions (phi, f, s);

  clear();
  draw_vof ("f", "s", edges = true, lw = 0.5);
  save ("vof.png");

  free(p);
}
