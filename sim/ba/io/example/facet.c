#include "grid/multigrid3D.h"
#include "fractions.h"
#include ".u/io/io.h"

#define LVL (6)

int
main(void)
{
  double r;
  origin(-0.5, -0.5, -0.5);
  init_grid(1 << LVL);
  scalar a[];
  r = 0.25;
  fraction(a, sq(x) + sq(y) + sq(z) - sq(r));
  DumpFacets(a, stdout);
}
