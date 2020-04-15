#include "grid/octree.h"
#include "utils.h"
#include "distance.h"
#include "fractions.h"
#include "view.h"

static const char *me = "distance";

int main()
{
  coord * p;
  coord min, max;
  double maxl;
  const char *cmd = "test -f distance.stl || "
      "(wget http://www.ai.sri.com/~reddy/book/models/horse.zip && "
      "unzip -o horse.zip && "
      "meshlabserver -i horse.ply -o distance.stl)";
  maxl = -HUGE;
  if (system (cmd) != 0) {
      fprintf(stderr, "%s: '%s' failed\n", me, cmd);
      exit(2);
  }
  p = input_stl (fopen ("distance.stl", "r"));
  bounding_box (p, &min, &max);
  foreach_dimension()
    if (max.x - min.x > maxl)
      maxl = max.x - min.x;
  init_grid (8);
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
  save ("horse.png");

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
}
