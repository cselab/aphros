changequote()dnl
function partstr_norm(i, j, u) {
    var nx, ny, n
    var X = 0, Y = 1
    p = Array(2)
    nx = (u[i+1][j+1]+2*u[i+1][j]+u[i+1][j-1]-u[i-1][j+1]-2*u[i-1][j]-u[i-1][j-1])/8
    ny = (u[i+1][j+1]-u[i+1][j-1]+2*u[i][j+1]-2*u[i][j-1]+u[i-1][j+1]-u[i-1][j-1])/8
    n =  -(Math.abs(nx) + Math.abs(ny))
    if (n == 0)
        throw new Error("n == 0")
    nx /= n
    ny /= n
    p[X] = nx
    p[Y] = ny
    return p
}

function _line(nx, ny, u) {
    u1 = 0.5 * nx / ny;
    if (u <= u1)
      return -0.5 * (nx + ny) + Math.sqrt(2*nx*ny*u)
    else
      return ny * (u - 0.5)
}

function partstr_line(nx, ny, u) {
    var t;
    if (ny < nx) {
        t = nx; nx = ny; ny = t
    }
    if (u < 0.5)
        return _line(nx, ny, u)
    else
        return -_line(nx, ny, 1 - u)
}

function partstr_vof_line(M, N, u, a) {
    const X = 0, Y = 1
    var n, u0, a0;
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++) {
            n = partstr_norm(i, j, u)
            u0 = u[i][j]
            a0 = partstr_line(n[X], n[Y], u0)
            a[i][j] = a0
        }
}
