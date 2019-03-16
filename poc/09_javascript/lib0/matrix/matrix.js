function matrix_read(file) {
    var a, nr, NR, l, line, i, j, x, m
    var lines = require('fs').readFileSync(file, 'ascii')
        .split('\n')
        .filter(Boolean)
    NR = lines.length
    if (NR < 0)
        throw new Error("no lines in file " + file)
    nr = 0
    line = lines[nr++]
    l = line.trim().split(" ")
    if (l.length != 2)
        throw new Error("expecting 'M N' got " + line)
    M = parseInt(l[0])
    if (Number.isNaN(M))
        throw new Error("expecting 'M', got " + l[0])
    N = parseInt(l[1])
    if (Number.isNaN(N))
        throw new Error("expecting 'N', got " + l[1])
    a = matrix_new(M, N)

    for (i = 0; i < M; i++) {
        if (nr >= NR)
            throw new Error("too few lines in ${file}, expecting ${NR}, got ${nr}")
        line = lines[nr++]
        l = line.trim().split(" ")
        if (l.length != N)
            throw new Error("expecting ${N} numbers, got '${line}'")
        for (j = 0; j < N; j++) {
            x = parseFloat(l[j])
            if (Number.isNaN(x))
                throw new Error("expecting float, got " + l[i])
            a[i][j] = x
        }
    }
    return a
}

function matrix_write(stream, M, N, a) {
    var i, j
    if (!Array.isArray(a))
        throw new Error("a is not an array: " + a)
    stream.write(M + " " + N + "\n")
    for (i = 0; i < M; i++) {
        if (a[i] === void 0)
            throw new Error("a[${i}] is undefined, M = ${M}, N = ${N}")
        for (j = 0; j < N; j++) {
            if (j > 0)
                stream.write(" ")
            stream.write(a[i][j] + "")
        }
        stream.write("\n")
    }
}

function matrix_lh_write(stream, m, n, M, N, a) {
    var i, j
    if (!Array.isArray(a))
        throw new Error("a is not an array: " + a)
    stream.write(M - m + " " + N - n + "\n")
    for (i = m; i < M; i++) {
        if (a[i] == void 0)
            throw new Error("a[${i}] is undefined, M = ${M}, N = ${N}")
        for (j = n; j < N; j++) {
            if (j > n) stream.write(" ")
            stream.write(a[i][j] + "")
        }
        stream.write("\n")
    }
}

function matrix_new(M, N) {
    var i, a
    a = new Array(M)
    for (i = 0; i < M; i++)
        a[i] = new Array(N)
    return a
}

function matrix_zero(M, N) {
    var i, j, a
    a = matrix_new(M, N)
    for (i = 0; i < M; i++)
        for (j = 0; j < M; j++)
            a[i][j] = 0
    return a
}

function matrix_halo_zero(M, N, h, a) {
    var i, j
    if (!Array.isArray(a))
        throw new Error("a is not an array: " + a)
    for (i = -h; i <     0; i++) a[i] = []
    for (i =  M; i < M + h; i++) a[i] = []

    for (i = -h; i < 0; i++)
        for (j = -h; j < N + h; j++) a[i][j] = 0;

    for (i = 0; i < M; i++) {
        for (j = -h; j < 0;      j++) a[i][j] = 0;
        for (j =   N; j < N + h; j++) a[i][j] = 0;
    }

    for (i = M; i < M + h; i++)
        for (j = -h; j < N + h; j++) a[i][j] = 0;
}

function matrix_transpose(M, N, a) {
    var i, j, b
    if (!Array.isArray(a))
        throw new Error("a is not an array: " + a)
    b = matrix_new(N, M)
    for (i = 0; i < M; i++)
        for (j = 0; j < M; j++)
            b[j][i] = a[i][j]
    return b
}
