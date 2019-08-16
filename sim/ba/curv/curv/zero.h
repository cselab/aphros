#include "fractions.h"
#include "curvature.h"

void curvature_zero(struct Curvature p)
{
  scalar c = p.c, kappa = p.kappa;
  foreach() {
    kappa[] = 0.;
  }
}

trace
cstats curvature_orig(struct Curvature p) {
  return curvature(p);
}
#define curvature curvature_zero
