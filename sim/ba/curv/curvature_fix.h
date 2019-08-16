#include "fractions.h"
#include "curvature.h"

trace
cstats curvature_fix (struct Curvature p)
{
  scalar c = p.c, kappa = p.kappa;
  double sigma = p.sigma ? p.sigma : 1.;
  int sh = 0, sf = 0, sa = 0, sc = 0;
  vector ch = c.height, h = automatic (ch);
  if (!ch.x.i)
    heights (c, h);

  /**
  On trees we set the prolongation and restriction functions for
  the curvature. */
  
#if TREE
  kappa.refine = kappa.prolongation = curvature_prolongation;
  kappa.restriction = curvature_restriction;
#endif

  /**
  We first compute a temporary curvature *k*: a "clone" of
  $\kappa$. */
  
  scalar k[];
  scalar_clone (k, kappa);

  foreach(reduction(+:sh) reduction(+:sf)) {

    /**
    If we are not in an interfacial cell, we set $\kappa$ to *nodata*. */

    if (!interfacial (point, c))
      k[] = nodata;

    /**
    Otherwise we try the standard HF curvature calculation first, and
    the "mixed heights" HF curvature second. */ 
    
    else if ((k[] = height_curvature (point, c, h)) != nodata)
      sh++;
    else if ((k[] = height_curvature_fit (point, c, h)) != nodata)
      sf++;
  }
  boundary ({k});

  foreach (reduction(+:sa) reduction(+:sc)) {
    
    /**
    We then construct the final curvature field using either the
    computed temporary curvature... */

    double kf;
    if (k[] < nodata)
      kf = k[];
    else if (interfacial (point, c)) {

      /**
      ...or the average of the curvatures in the $3^{d}$ neighborhood
      of interfacial cells. */
      
      double sk = 0., a = 0.;
      foreach_neighbor(1)
	if (k[] < nodata)
	  sk += k[], a++;
      if (a > 0.)
	kf = sk/a, sa++;
      else

	/**
	Empty neighborhood: we try centroids as a last resort. */

	kf = centroids_curvature_fit (point, c), sc++;
    }
    else
      kf = nodata;

    /**
    We add or set *kappa*. */
    
    if (kf == nodata)
      kappa[] = nodata;
    else if (p.add)
      kappa[] += sigma*kf;
    else
      kappa[] = sigma*kf;      
  }
  boundary ({kappa});

  return (cstats){sh, sf, sa, sc};
}

trace
cstats curvature_orig(struct Curvature p) {
  return curvature(p);
}
#define curvature curvature_fix
