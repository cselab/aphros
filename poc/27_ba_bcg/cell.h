/**
# Bell-Collela-Glaz advection scheme

The function below implements the 2nd-order, unsplit, upwind scheme of
[Bell-Collela-Glaz, 1989](references.bib#bell89). Given a centered
scalar field *f*, a face vector field *uf* (possibly weighted by a
face metric), a timestep *dt* and a source term field *src*, it fills
the face vector field *flux* with the components of the advection
fluxes of *f*. */

void tracer_fluxes (scalar f,
    face vector uf,
    face vector flux,
    double dt,
    (const) scalar src)
{

  /**
    We first compute the cell-centered gradient of *f* in a locally-allocated
    vector field. */

  vector g[];
  gradients ({f}, {g});

  /**
    For each face, the flux is composed of two parts... */

  foreach_face() {

    /**
      A normal component... (Note that we cheat a bit here, `un` should
      strictly be `dt*(uf.x[i] + uf.x[i+1])/((fm.x[] +
      fm.x[i+1])*Delta)` but this causes trouble with boundary
      conditions (when using narrow '1 ghost cell' stencils)). */

    double un = dt*uf.x[]/(fm.x[]*Delta + SEPS), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = f[i] + (src[] + src[-1])*dt/4. + s*(1. - s*un)*g.x[i]*Delta/2.;

    /**
      and tangential components... */

    double vn = (uf.y[i] + uf.y[i,1])/(fm.y[i] + fm.y[i,1]);
    double fyy = vn < 0. ? f[i,1] - f[i] : f[i] - f[i,-1];
    f2 -= dt*vn*fyy/(2.*Delta);

    double wn = (uf.z[i] + uf.z[i,0,1])/(fm.z[i] + fm.z[i,0,1]);
    double fzz = wn < 0. ? f[i,0,1] - f[i] : f[i] - f[i,0,-1];
    f2 -= dt*wn*fzz/(2.*Delta);


    flux.x[] = f2*uf.x[];
  }

  /**
    Boundary conditions ensure the consistency of fluxes across
    variable-resolution boundaries (on adaptive meshes). */

  boundary_flux ({flux});
}

/**
  The function below uses the *tracer_fluxes* function to integrate the
  advection equation, using an explicit scheme with timestep *dt*, for
  each tracer in the list. */

struct Advection {
  scalar * tracers;
  face vector u;
  double dt;
  scalar * src; // optional
};

void advection (struct Advection p)
{

  /**
    If *src* is not provided we set all the source terms to zero. */

  scalar * lsrc = p.src;
  if (!lsrc) {
    const scalar zero[] = 0.;
    for (scalar s in p.tracers)
      lsrc = list_append (lsrc, zero);
  }

  assert (list_len(p.tracers) == list_len(lsrc));
  scalar f, src;
  for (f,src in p.tracers,lsrc) {
    face vector flux[];
    tracer_fluxes (f, p.u, flux, p.dt, src);
    foreach()
      foreach_dimension()
      f[] += p.dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
  }
  boundary (p.tracers);

  if (!p.src)
    free (lsrc);
}
