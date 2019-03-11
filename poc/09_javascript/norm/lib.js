function InitGrid(nx, ny, u) {
  var b = 2
  for (var j = -b; j < ny + b; j++) {
      u[j] = [];
      for (var i = -b; i < nx + b; i++) {
          u[j][i] = (i * j) / (nx * ny)
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

function Norm(i, j, u, /**/ p) {
    var nx, ny, n
    var X = 0, Y = 1
    nx = (u[i+1][j+1]+2*u[i+1][j]+u[i+1][j-1]-u[i-1][j+1]-2*u[i-1][j]-u[i-1][j-1])/8
    ny = (u[i+1][j+1]-u[i+1][j-1]+2*u[i][j+1]-2*u[i][j-1]+u[i-1][j+1]-u[i-1][j-1])/8
    n =  -(Math.abs(nx) + Math.abs(ny))
    nx /= n
    ny /= n
    p[X] = nx
    p[Y] = ny
}
