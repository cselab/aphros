/* static */
function _sq(a) { return a*a }
function _E(a) { return [Math.cos(a), Math.sin(a)] }
function _axpy(a, x, y,   /**/ b) {
    const X = 0, Y = 1
    b[X] = a*x[X] + y[X]
    b[Y] = a*x[Y] + y[Y]
}
function _line(nx, ny, u) {
    var u1
    u1 = 0.5 * nx / ny
    if (u <= u1)
        return -0.5 * (nx + ny) + Math.sqrt(2*nx*ny*u)
    else
        return ny * (u - 0.5)
}
function _msg(s) { return process.stderr.write(s) }
/**/

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

function partstr_ends_read(file) {
    return matrix_read(file)
}

function partstr_ends_write(stream, ends) {
    var n
    n = ends.length
    matrix_write(stream, n, 4, ends)
}

function partstr_cell_ends_gnuplot_write(stream, n, ends) {
    const AX = 0, AY = 1, BX = 2, BY = 3
    const X = 0, Y = 1
    var i, e
    for (i = 0; i < n; i++) {
        e = ends[i]
        if (i > 0)
            stream.write("\n")
        stream.write(`${e[AX]} ${e[AY]}\n`)
        stream.write(`${e[BX]} ${e[BY]}\n`)
    }
}

function partstr_ends_gnuplot_write(stream, M, N, ends) {
    const AX = 0, AY = 1, BX = 2, BY = 3
    const X = 0, Y = 1
    var i, e
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++) {
            e = ends[i][j]
            if (i != 0 || j != 0)
                stream.write("\n")
            stream.write(`${e[AX]} ${e[AY]}\n`)
            stream.write(`${e[BX]} ${e[BY]}\n`)
        }
}

function partstr_part(nh, hp, p, a, t, /**/ xx) {
    var n, j, jp
    if (!Array.isArray(xx))
        throw new Error(`xx is not an array: ${xx}`)
    n = 2*nh + 1
    xx[nh] = p.slice()
    for (j = 0; j < nh; j++) {
        jp = j + 0.5
        _axpy( hp, _E(a + t*jp), xx[nh + j], /**/ xx[nh + j + 1])
        _axpy(-hp, _E(a - t*jp), xx[nh - j], /**/ xx[nh - j - 1])
    }
}

function partstr_segcirc(k, l, d) {
    var t1, t2
    t1 = 1 - _sq(k*l)
    t2 = 1 - _sq(k*d)
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

    mdc = _sq(u) + _sq(v)
    mdc = Math.sqrt(mdc)/2.0

    dc =  _sq(p) + _sq(q)
    dc = Math.sqrt(dc)

    g = _sq(u) + _sq(v)
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
    g = _sq(u) + _sq(v)
    if (g > 0)
        k /= g

    if (k < 0) k = 0
    if (k > 1) k = 1

    return [a[X] + k*u, a[Y] + k*v]
}

function partstr_nearest_ends(n, ends, x, k) {
    const X = 0, Y = 1
    const AX = 0, AY = 1, BX = 2, BY = 3
    var i, a, b, e, y, z, d
    var m, j
    if (!Array.isArray(ends))
        throw new Error(`ends is not an array: ${ends}`)
    a = Array(2); b = Array(2)
    for (i = 0; i < n; i++) {
        e = ends[i]
        if (!Array.isArray(e))
            throw new Error(`e is not an array: ${e}`)
        a[X] = e[AX]; a[Y] = e[AY]
        b[X] = e[BX]; b[Y] = e[BY]
        y = partstr_nearest(a, b, x)
        d = _sq(x[X] - y[X]) + _sq(x[Y] - y[Y])
        if (i == 0 || d < m) {
            m = d
            j = i
            z = partstr_shsegcirc(k, a, b, y)
        }
    }
    return z
}

function partstr_force(ne, ends, np, xx, k, eta, /**/ ff) {
    const X = 0, Y = 1
    var i, x, y
    if (!Array.isArray(ends))
        throw new Error(`ends is not an array: ${ends}`)
    if (!Array.isArray(xx))
        throw new Error(`xx is not an array: ${xx}`)
    if (!Array.isArray(ff))
        throw new Error(`ff is not an array: ${ff}`)
    for (i = 0; i < np; i++) {
        x = xx[i]
        y = partstr_nearest_ends(ne, ends, x, k)
        ff[i][X] = eta*(y[X] - x[X])
        ff[i][Y] = eta*(y[Y] - x[Y])
    }
}

function partstr_force_write(stream, n, xx, ff) {
    const X = 0, Y = 1
    var i, x, f
    if (!Array.isArray(ff))
        throw new Error(`ff is not an array: ${ff}`)
    if (!Array.isArray(xx))
        throw new Error(`xx is not an array: ${xx}`)

    for (i = 0; i < n; i++) {
        x = xx[i]; f = ff[i]
        stream.write(`${x[X]} ${x[Y]} ${f[X]} ${f[Y]}\n`)
    }
}

function partstr_dxda(nh, hp, a, t, /**/ xx) {
    var n, j, jp, pi
    if (!Array.isArray(xx))
        throw new Error(`xx is not an array: ${xx}`)
    n = 2*nh + 1
    pi = Math.PI
    xx[nh] = [0, 0]
    for (j = 0; j < nh; j++) {
        jp = j + 0.5
        _axpy( hp, _E(a + t*jp + pi/2), xx[nh + j], /**/ xx[nh + j + 1])
        _axpy(-hp, _E(a - t*jp + pi/2), xx[nh - j], /**/ xx[nh - j - 1])
    }
}

function partstr_dxdt(nh, hp, a, t, /**/ xx) {
    var n, j, jp, pi
    if (!Array.isArray(xx))
        throw new Error(`xx is not an array: ${xx}`)
    n = 2*nh + 1
    pi = Math.PI
    xx[nh] = [0, 0]
    for (j = 0; j < nh; j++) {
        jp = j + 0.5
        _axpy(hp*jp, _E(a + t*jp + pi/2), xx[nh + j], /**/ xx[nh + j + 1])
        _axpy(hp*jp, _E(a - t*jp + pi/2), xx[nh - j], /**/ xx[nh - j - 1])
    }
}

function partstr_curv(hp, t) {
    return Math.sqrt(2)/hp*Math.sin(t)/Math.sqrt(1.0 + Math.cos(t))
}

function _minus(n, a, b, /**/ c) {
    const X = 0, Y = 1
    var i
    for (i = 0; i < n; i++) {
        c[i][X] = a[i][X] - b[i][X]
        c[i][Y] = a[i][Y] - b[i][Y]
    }
}
function _substr(n, a, /**/ b) {
    const X = 0, Y = 1
    var i
    for (i = 0; i < n; i++) {
        b[i][X] -= a[i][X]
        b[i][Y] -= a[i][Y]
    }
}
function _append(n, a, b /**/) {
    const X = 0, Y = 1
    var i
    for (i = 0; i < n; i++) {
        b[i][X] += a[i][X]
        b[i][Y] += a[i][Y]
    }
}
function _dot(n, a, b) {
    const X = 0, Y = 1
    var i, d
    d = 0
    for (i = 0; i < n; i++) {
        d += a[i][X]*b[i][X]
        d += a[i][Y]*b[i][Y]
    }
    return d
}
function partstr_step(nh, ff, eta, hp, /*io*/ State) {
    const X = 0, Y = 1
    var p, a, t, x0, x1, f, n, dx, dd

    p = State.p.slice()
    a = State.a
    t = State.t

    f = ff[nh]
    n = 2*nh + 1

    x0 = matrix_new(n, 2) /**/
    x1 = matrix_new(n, 2) /* todo */
    dx = matrix_new(n, 2) /**/
    dd = matrix_new(n, 2) /**/

    /* p */
    partstr_part(nh, hp, p, a, t, /**/ x0)
    p[X] += f[X]
    p[Y] += f[Y]
    partstr_part(nh, hp, p, a, t, /**/ x1)
    _minus(n, x1, x0, /**/ dx)
    _append(n, dx, /**/ x0)

    /* a */
    partstr_dxda(nh, hp, a, t, /**/ dd)
    a += _dot(n, ff, dd)/_dot(n, dd, dd)
    partstr_part(nh, hp, p, a, t, /**/ x1)
    _minus(n, x1, x0, /**/ dx)

    /* t */
    partstr_dxdt(nh, hp, a, t, /**/ dd)
    t += _dot(n, ff, dd)/_dot(n, dd, dd)

    State.p = p.slice()
    State.a = a
    State.t = t
}

function Partstr(nh, hp, eta) {
    const n = 2*nh + 1
    this.nh = nh
    this.hp = hp
    this.eta = eta
    this.ff = matrix_new(n, 2)
    this.xx = matrix_new(n, 2)

    this.start = function(ne, ends, a, t, p) {
        if (!Array.isArray(ends))
            throw new Error(`ends is not an array: ${ends}`)
        if (!Array.isArray(p))
            throw new Error(`p is not an array: ${p}`)
        this.ne = ne
        this.ends = ends
        this.a = a
        this.t = t
        this.p = p.slice()
    }

    this.step = function() {
        var nh, p, a, t, hp, eta, ends, ff, State, k, ne, xx
        nh = this.nh
        const n = 2*nh + 1
        hp = this.hp
        eta = this.eta
        ends = this.ends
        ne = this.ne
        xx = this.xx
        ff = this.ff
        p = this.p; a = this.a; t = this.t
        k = partstr_curv(hp, t)
        partstr_part(nh, hp, p, a, t, /**/ xx)
        partstr_force(ne, ends, n, xx, k, eta,     ff)
        State = { p: p, a: a, t: t }
        partstr_step(nh, ff, eta, hp, /*io*/ State)
        this.a = State.a; this.t = State.t; this.p = State.p
        this.k = k
        this.xx = xx
        this.ff = ff
    }
}
