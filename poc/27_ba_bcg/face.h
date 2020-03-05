void prediction()
{
  vector du;
  foreach_dimension() {
    scalar s = new scalar;
    du.x = s;
  }

  // du : approximation of slope

  if (u.x.gradient)
    foreach()
      foreach_dimension() {
        du.x[] = u.x.gradient (u.x[-1], u.x[], u.x[1])/Delta;
      }
  else
    foreach()
      foreach_dimension() {
        du.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
      }
  boundary ((scalar *){du});

  trash ({uf});
  foreach_face() {
    double un = dt*(u.x[] + u.x[-1])/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    uf.x[] = u.x[i] + (g.x[] + g.x[-1])*dt/4. + s*(1. - s*un)*du.x[i]*Delta/2.;

    double fyy = u.y[i] < 0. ? u.x[i,1] - u.x[i] : u.x[i] - u.x[i,-1];
    uf.x[] -= dt*u.y[i]*fyy/(2.*Delta);

    double fzz = u.z[i] < 0. ? u.x[i,0,1] - u.x[i] : u.x[i] - u.x[i,0,-1];
    uf.x[] -= dt*u.z[i]*fzz/(2.*Delta);

    uf.x[] *= fm.x[];
  }
  boundary ((scalar *){uf});

  delete ((scalar *){du});
}


/*


Linear profile:

  u_x = du / dx


We first define an initial approximation

  du = (u[i+1] - u[i-1]) / 2

or at boundaries

  du =  (u[2] - u[1] - 2 * u[0.5]) / 2

which will then be limited.

Limiter:

|du| <= 2 * max(|u[i+1] - u[i]|, |u[i] - u[i-1]|)

and du=0 if (u[i+1] - u[i]) * (u[i] - u[i-1]) < 0

or at boundaries such that the linear profile evaluated at the boundary lies
between the cell average and the specified boundary value.

   */



// prediction() if g=0 and un>0

void prediction_g0() {
  double un = dt * (ux[] + ux[-1]) / (2. * Delta);
  int i = -1;
  uf.x[] = ux[i] + (1. - un) * dux[i] * Delta / 2.;

  double fyy = uy[i] < 0. ? ux[i, 1] - ux[i] : ux[i] - ux[i, -1];
  uf.x[] -= dt * uy[i] * fyy / (2. * Delta);

  double fzz = uz[i] < 0. ? ux[i, 0, 1] - ux[i] : ux[i] - ux[i, 0, -1];
  uf.x[] -= dt * uz[i] * fzz / (2. * Delta);

  uf.x[] *= fm.x[];
}

void cell() {
  g = gradient(f);

  double s = sign(dt * uf.x[]);
  int i = -(s + 1.) / 2.;

  // cell value of field
  double f = ff_f[i];

  // advance with spatial derivative in normal direction
  f += s * g.x[i] * Delta / 2.;

  // temporal derivative, source
  double ut = (src[] + src[-1]) / 2;

  // convective terms
  double un = uf.x[] / fm.x[];
  double gx = g.x[i];
  ut -= un * gx;

  double vn = (uf.y[i] + uf.y[i, 1]) / (fm.y[i] + fm.y[i, 1]);
  double gy = (vn < 0. ? f[i, 1] - f[i] : f[i] - f[i, -1]) / Delta;
  ut -= vn * gy;

  double wn = (uf.z[i] + uf.z[i, 0, 1]) / (fm.z[i] + fm.z[i, 0, 1]);
  double gz = (wn < 0. ? f[i, 0, 1] - f[i] : f[i] - f[i, 0, -1]) / Delta;
  ut -= wn * gz;

  // advance with temporal derivative
  f += ut * dt / 2;


  // Compared to BCG 1989,
  // this computation does not include limiters,
  // and upwinding based on the Riemann problem.
}

