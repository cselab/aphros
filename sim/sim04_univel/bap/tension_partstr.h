/**
# Surface tension

Surface tension can be expressed as the interfacial force density
$$
\phi\nabla f
$$
with $f$ the volume fraction describing the interface and the potential
$$
\phi = \sigma\kappa
$$
with $\sigma$ the (constant) surface tension coefficient and $\kappa$
the interface mean curvature. */

#include "iforce.h"
#include "partstr.h"

/**
The surface tension coefficient is associated to each VOF tracer. */

attribute {
  double sigma;
}

/**
## Stability condition

The surface tension scheme is time-explicit so the maximum timestep is
the oscillation period of the smallest capillary wave.
$$
T = \sqrt{\frac{\rho_{m}\Delta_{min}^3}{\pi\sigma}}
$$
with $\rho_m=(\rho_1+\rho_2)/2.$ and $\rho_1$, $\rho_2$ the densities
on either side of the interface. */

event stability (i++) {

  /**
  We first compute the minimum and maximum values of $\alpha/f_m =
  1/\rho$, as well as $\Delta_{min}$. */

  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin)) {
    if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
    if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
    if (Delta < dmin) dmin = Delta;
  }
  double rhom = (1./amin + 1./amax)/2.;

  /**
  The maximum timestep is set using the sum of surface tension
  coefficients. */

  double sigma = 0.;
  for (scalar c in interfaces)
    sigma += c.sigma;
  if (sigma) {
    double dt = sqrt (rhom*cube(dmin)/(pi*sigma));
    if (dt < dtmax)
      dtmax = dt;
  }
}

/**
## Definition of the potential

We overload the acceleration event to define the potential
$\phi=\sigma\kappa$. */

event acceleration (i++)
{
  
  /**
  We check for all VOF interfaces for which $\sigma$ is non-zero. */

  for (scalar f in interfaces)
    if (f.sigma) {
      
      /**
      If $\phi$ is already allocated, we add $\sigma\kappa$, otherwise
      we allocate a new field and set it to $\sigma\kappa$. */

      scalar phi = f.phi;
      if (phi.i)
	curvature_partstr(f, phi, f.sigma, add = true);
      else {
	phi = new scalar;
	curvature_partstr(f, phi, f.sigma, add = false);
	f.phi = phi;
      }
    }
}
