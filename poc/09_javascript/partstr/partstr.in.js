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
    var u1
    u1 = 0.5 * nx / ny
    if (u <= u1)
      return -0.5 * (nx + ny) + Math.sqrt(2*nx*ny*u)
    else
      return ny * (u - 0.5)
}

function partstr_line(nx, ny, u) {
    var t
    nx = Math.abs(nx)
    ny = Math.abs(ny)
    if (ny < nx) {
        t = nx; nx = ny; ny = t
    }
    if (u < 0.5)
        return _line(nx, ny, u)
    else
        return -_line(nx, ny, 1 - u)
}

function partstr_vof_line(M, N, u, /**/ a) {
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

function partstr_vof_ends(M, N, u, /**/ ends) {
    const X = 0, Y = 1
    const AX = 0, AY = 1, BX = 2, BY = 3
    const h = 1.0
    var n, a, i, j;
    var e = Array(4)
    for (i = 0; i < M; i++)
    for (j = 0; j < N; j++) {
        n = partstr_norm(i, j, u)
        a = partstr_line(n[X], n[Y], u[i][j])
        partstr_ends(n[X], n[Y], a, /**/ e)
        e[AX] += i*h + h/2; e[AY] += j*h + h/2
        e[BX] += i*h + h/2; e[BY] += j*h + h/2
        ends[i][j] = e.slice()
    }
}

function partstr_ends(nx, ny, a, /**/ e) {
    var x, y, u, v, h
    var h  = 0.5

    x = (a + h*ny)/nx
    y = (a + h*nx)/ny

    u = (a - h*ny)/nx
    v = (a - h*nx)/ny

    e[0] = e[1] = e[2] = e[3] = 0
    j = 0
    if (-h <= x && x <= h) {
        e[2*j] = x; e[2*j+1] = -h; j++
    }

    if (-h <= u && u <= h) {
        e[2*j] = u; e[2*j + 1] = h; j++
    }

    if (j < 2 && -h <= y && y <= h) {
        e[2*j] = -h; e[2*j + 1] = y; j++
    }

    if (j < 2 && -h <= v && v <= h) {
        e[2*j] = h; e[2*j + 1] = v; j++
    }

    if (j == 1) {
        e[2*j] = e[0]; e[2*j + 1] = e[1]
    }
}
