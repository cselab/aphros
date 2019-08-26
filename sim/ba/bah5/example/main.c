#include <h5.h>
#include "grid/multigrid3D.h"
#include ".u/bah5.h"

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
  bah5_list({a, b}, "%s_000"); /* '%s is replaced by a field name */
  bah5_list0(scalars, names, "%s_000");
}
