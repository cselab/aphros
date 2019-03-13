function matrix_read(file) {
    var a, nr, NR, l, line, i, j, x, m
    var lines = require('fs').readFileSync(file, 'ascii')
        .split('\n')
        .filter(Boolean)
    NR = lines.length
    if (NR < 0)
        throw new Error(`no lines in file ${file}`)
    nr = 0
    line = lines[nr++]
    l = line.trim().split(" ")
    if (l.length != 2)
        throw new Error(`expecting 'M N' got '${line}'`)
    M = parseInt(l[0])
    if (Number.isNaN(M))
        throw new Error(`expecting 'M', got '${l[0]}'`)
    N = parseInt(l[1])
    if (Number.isNaN(N))
        throw new Error(`expecting 'N', got '${l[1]}'`)
    a = matrix_new(M, N)

    for (i = 0; i < M; i++) {
        if (nr >= NR)
            throw new Error(`too few lines in ${file}, expecting ${NR}, got ${nr}`)
        line = lines[nr++]
        l = line.trim().split(" ")
        if (l.length != N)
            throw new Error(`expecting ${N} numbers, got '${line}'`)
        for (j = 0; j < N; j++) {
            x = parseFloat(l[j])
            if (Number.isNaN(x))
                throw new Error(`expecting float, got '${l[i]}'`)
            a[i][j] = x
        }
    }
    return a
}

function matrix_write(stream, M, N, a) {
    var i, j
    stream.write(`${M} ${N}\n`)
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            if (j > 0)
                stream.write(" ")
            stream.write(`${a[i][j]}`)
        }
        stream.write("\n")
    }
}

function matrix_lh_write(stream, m, n, M, N, a) {
    var i, j
    for (i = m; i < M; i++) {
        if (a[i] == void 0)
            throw new Error(`a[${i}] is undefined`)
        for (j = n; j < N; j++) {
            if (j > n) stream.write(" ")
            stream.write(`${a[i][j]}`)
        }
        stream.write("\n")
    }
}

function matrix_new(M, N) {
    var i
    a = new Array(M)
    for (i = 0; i < M; i++)
        a[i] = new Array(N)
    return a
}

function matrix_zero(M, N) {
    var i, j
    a = matrix_new(M, N)
    for (i = 0; i < M; i++)
        for (j = 0; j < M; j++)
            a[i][j] = 0
    return a
}

function matrix_halo_zero(M, N, h, a) {
    var i, j

    for (i = -h; i <     0; i++) a[i] = []
    for (i =  N; i < N + h; i++) a[i] = []

    for (i = -h; i < 0; i++)
        for (j = -h; j < N + h; j++) a[i][j] = 0;

    for (i = 0; i < M; i++) {
        for (j = -h; j < 0;      j++) a[i][j] = 0;
        for (j =   N; j < N + h; j++) a[i][j] = 0;
    }

    for (i = M; i < M + h; i++)
        for (j = -h; j < N + h; j++) a[i][j] = 0;
}
function partstr_norm(i, j, u) {
    var nx, ny, n
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
    var x, y, u, v, h, cross, t
    var h  = 0.5

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
    var e, ans, k
    ans = Array(); k = 1
    for (i = m - s; i < m + s + 1; i++) {
	if (i <  0) continue
	if (i >= M) continue
	for (j = n - s; j < n + s + 1; j++) {
	    if (j <  0) continue
	    if (j >= N) continue
	    e = ends[i][j]
	    if (i == n && j == m)
		ans[0] = e
	    else
		ans[k++] = e
	}
    }
    return ans
}

function partstr_ends_write(stream, ends) {
    var n
    n = e.length
    matrix_write(stream, n, 4, ends)
}


const AX = 0, AY = 1, BX = 2, BY = 3
const stdout = process.stdout
const stderr = process.stdout
const hl = 2

var i = 2, j, e
var file = process.argv[i++]
var u = matrix_read(file)
var M = u.length
var N = u[0].length

matrix_halo_zero(M, N, hl, u)

var ends = matrix_new(M, N)
partstr_vof_ends(M, N, u, ends)

i = 2, j = 2
e = partstr_cell_ends(M, N, i, j, ends)

/*
n = e.length
for (i = 0; i < 1; i++) {
    f = e[i]
    stdout.write(`${f[AX]} ${f[AY]}\n`)
    stdout.write(`${f[BX]} ${f[BY]}\n`)
    stdout.write("\n")
    }*/

partstr_ends_write(stdout, e)
