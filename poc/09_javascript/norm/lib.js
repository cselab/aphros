function InitGrid(nx, ny, u) {
  for (var j = 0; j < ny; j++) {
      u[j] = [];
      for (var i = 0; i < nx; i++) {
	  u[j][i] = 1. / (sqr(i - (nx - 1) * 0.5) +
		    sqr(j - (ny - 1) * 0.5));
	  if (u[j][i] < 0.3) {
	      u[j][i] = 0.;
	  }
	  u[j][i] = Clip(u[j][i], 0., 1.)
      }
  }
  return u;
}

function Clip(a,l,h) {
    return Math.min(Math.max(a, l), h);
}

function sqr(a) {
    return a * a;
}

