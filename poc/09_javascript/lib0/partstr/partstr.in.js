changequote()dnl
function partstr_norm(i, j, u) {
    var nx, ny, n, p
    var X = 0, Y = 1
    p = Array(2)
    nx = (u[i+1][j+1]+2*u[i+1][j]+u[i+1][j-1]-u[i-1][j+1]-2*u[i-1][j]-u[i-1][j-1])/8
    ny = (u[i+1][j+1]-u[i+1][j-1]+2*u[i][j+1]-2*u[i][j-1]+u[i-1][j+1]-u[i-1][j-1])/8
    n =  -(Math.abs(nx) + Math.abs(ny) + 1e-10)
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
    var i, j
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

function partstr_ends(M, N, a, /**/ e) {
    const AX = 0, AY = 1, BX = 2, BY = 3
    const h  = 0.5
    var x, y, u, v, cross, t
    var j

    x = (a + h*N)/M
    y = (a + h*M)/N

    u = (a - h*N)/M
    v = (a - h*M)/N

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

    cross = (e[BY]-e[AY])*M-(e[BX]-e[AX])*N
    if (cross < 0) {
        t = e[AX]; e[AX] = e[BX]; e[BX] = t
        t = e[AY]; e[AY] = e[BY]; e[BY] = t
    }
}

function partstr_cell_ends(M, N, m, n, ends) {
    const s = 1
    var e, ans, i, j, k, Seen
    Seen = false
    ans = Array(); k = 1
    for (i = m - s; i < m + s + 1; i++) {
        if (i <  0) continue
        if (i >= M) continue
        for (j = n - s; j < n + s + 1; j++) {
            if (j <  0) continue
            if (j >= N) continue
            e = ends[i][j]
            if (i == m && j == n) {
                Seen = true
                ans[0] = e
            } else
                ans[k++] = e
        }
    }
    if (!Seen)
        throw new Error(
            `have not seen ends[${m}][${n}], M = ${M}, N = ${N}\n`)
    return ans
}

function partstr_ends_write(stream, ends) {
    var n
    n = e.length
    matrix_write(stream, n, 4, ends)
}

function partstr_ends_read(stream) {
    return matrix_read(stream, n, 4, ends)
}

function _E(a) { return [Math.cos(a), Math.sin(a)] }
function _axpy(a, x, y,   /**/ b) {
    const X = 0, Y = 1
    b[X] = a*x[X] + y[X]
    b[Y] = a*x[Y] + y[Y]
}
function partstr_part(nh, hp, p, a, t) {
    var n, j, jp, xx
    n = 2*nh + 1
    xx = matrix_new(n, 2)
    xx[nh] = p.slice()
    for (j = 0; j < nh; j++) {
        jp = j + 0.5
        _axpy( hp, _E(a + t*jp), xx[nh + j], /**/ xx[nh + j + 1])
        _axpy(-hp, _E(a - t*jp), xx[nh - j], /**/ xx[nh - j - 1])
    }
    return xx
}

function partstr_segcirc(k, l, d) {
    var t1, t2
    t1 = 1 - (k*l)**2
    t2 = 1 - (k*d)**2
    if (t1 < 0)
        throw new Error(`t1=${t1} < 0, k=${k}, l=${l}, d=${d}\n`)
    if (t2 < 0)
        throw new Error(`t2=${t2} < 0, k=${k}, l=${l}, d=${d}\n`)
    t1 = Math.sqrt(t1)
    t2 = Math.sqrt(t2)
    if (t1 + t2 === 0.0)
        throw new Error(`t1 + t2 === 0, k=${k}, l=${l}, d=${d}\n`)
    return k*(l*l - d*d)/(t1 + t2)
}

function partstr_shsegcirc(k, a, b, x) {
    const X = 0, Y = 1
    const sc = 1
    var u, v, p, q, g, dc, mdc, s
    u = b[X] - a[X]
    v = b[Y] - a[Y]
    p = 0.5*(b[X] + a[X]) - x[X]
    q = 0.5*(b[Y] + a[Y]) - x[Y]

    mdc = u**2 + v**2
    mdc = Math.sqrt(mdc)/2.0

    dc =  p**2 + q**2
    dc = Math.sqrt(dc)

    g = u**2 + v**2
    if (g > 0) {
        g = Math.sqrt(g)
        u /= g
        v /= g
    }

    s = sc*partstr_segcirc(k, mdc, dc);
    return [x[X] + s*v, x[Y] - s*u]
}

function partstr_nearest(a, b, x) {
    const X = 0, Y = 1
    var u, v, p, q, k, g

    u = b[X] - a[X]
    v = b[Y] - a[Y]

    p = x[X] - a[X]
    q = x[Y] - a[Y]

    k = u*p + v*q
    g = u**2 + v**2
    if (g > 0)
        k /= g

    if (k < 0) k = 0
    if (k > 1) k = 1

    return [a[X] + k*u, a[Y] + k*v]
}

function partstr_nearest_ends(n, ends, x, k) {
    const X = 0, Y = 1
    const AX = 0, AY = 1, BX = 2, BY = 3
    var i, a, b, e, y, d
    var m, j
    a = Array(2); b = Array(2)
    for (i = 0; i < n; i++) {
        e = ends[i]
        a[X] = e[AX]; a[Y] = e[AY]
        b[X] = e[BX]; b[Y] = e[BY]
        y = partstr_nearest(a, b, x)
        d = (x[X] - y[X])**2 + (x[Y] - y[Y])**2
        if (i == 0 || d > m) {
            j = i
            m = d
        }
    }
    y = partstr_shsegcirc(k, a, b, y)
    return y
}
