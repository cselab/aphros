#include "grid/multigrid3D.h"
#include ".u/io/io.h"

#define LVL (6)

int
main(void)
{
  init_grid(1 << LVL);
  scalar a[], b[];
  char *names[] = {"pressure", "density"};
  scalar *scalars = {a, b};
  foreach() {
    a[] = (x - X0) * z;
    b[] = (x - X0);
  }
  io_write({a, b}, stdout);
}
