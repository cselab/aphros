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

function partstr_cell_ends(M, N, m, n, /**/ ends) {
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
    var n, j, jp
    n = 2*nh + 1
    xx = matrix_new(n, 2)
    xx[nh] = p.slice()
    for (j = 0; j < nh; j++) {
	jp = j + 0.5
	process.stderr.write(`${a} ${t} ${jp} ${_E(a + t*jp)}\n`)
	_axpy(hp, _E(a + t*jp), xx[nh + j], /**/ xx[nh + j + 1])
	_axpy(hp, _E(a - t*jp), xx[nh - j], /**/ xx[nh - j - 1])
    }
    process.stderr.write(`${xx}\n`)
}
