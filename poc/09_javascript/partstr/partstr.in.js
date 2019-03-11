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
    u1 = 0.5 * nx / ny
    if (u <= u1)
      return -0.5 * (nx + ny) + Math.sqrt(2*nx*ny*u)
    else
      return ny * (u - 0.5)
}

function partstr_line(nx, ny, u) {
    var t
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
    var n, u0, a0
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++) {
            n = partstr_norm(i, j, u)
            u0 = u[i][j]
            a0 = partstr_line(n[X], n[Y], u0)
            a[i][j] = a0
        }
}

function partstr_ends(nx, ny, a, /**/ e) {
    const X = 0, Y = 1
    const AX = 0, AY = 1, BX = 2, BY = 3
    var xl = Array(2)
    var xr = Array(2)
    var  n = Array(2)
    var hh = Array(2)

    n[X] = nx
    n[Y] = ny

    xl[X] = (a + 0.5 * n[1]) / n[0]
    xl[Y] = (a + 0.5 * n[0]) / n[1]

    xr[X] = (a - 0.5 * n[1]) / n[0]
    xr[Y] = (a - 0.5 * n[0]) / n[1]

    e[AX] = e[AY] = e[BX] = e[BY] = 0

    hh[X] = hh[Y] = 0.5

    j = 0
    if (-hh[0] <= xl[0] && xl[0] <= hh[0]) {
	e[2*j]   = xl[0]
	e[2*j+1] = -hh[1]
	j++
    }
    
    if (-hh[0] <= xr[0] && xr[0] <= hh[0]) {
	e[2*j] = xr[0]
	e[2*j + 1] = hh[1]
	j++
    }
    
    if (i < 2 && -hh[1] <= xl[1] && xl[1] <= hh[1]) {
	e[2*j] = -hh[0]
	e[2*j + 1] = xl[1]
	j++
    }
    
    if (i < 2 && -hh[1] <= xr[1] && xr[1] <= hh[1]) {
	e[2*j] = hh[0]
	e[2*j + 1] = xr[1]
	j++
    }
    
    if (j == 1) {
	e[2*j] = e[0]
	e[2*j + 1] = e[1]
    }
    return e;

    process.stderr.write(`xl: ${xl[X]} ${xl[Y]}\n`)
    process.stderr.write(`xr: ${xr[X]} ${xr[Y]}\n`)
}
