changequote()

function matrix_read(file) {
    var a, nr, NR, l, line, i, j, x, m
    var lines = require('fs').readFileSync(file, 'ascii')
        .split('\n')
        .filter(Boolean);

    NR = lines.length
    if (NR < 0)
        throw new Error(`no lines in file ${file}`)

    nr = 0
    line = lines[nr++]
    l = line.split(" ")
    if (l.length != 2)
        throw new Error(`expecting 'M N' got '${line}'`)
    M = parseInt(l[0])
    if (Number.isNaN(M))
        throw new Error(`expecting 'M', got '${l[0]}'`)
    N = parseInt(l[1])
    if (Number.isNaN(N))
        throw new Error(`expecting 'N', got '${l[1]}'`)
    a = new Array(M)
    for (i = 0; i < M; i++)
        a[i] = new Array(N)

    for (i = 0; i < M; i++) {
        if (nr >= NR)
            throw new Error(`too few lines in ${file}, expecting ${NR}, got ${nr}`)
        line = lines[nr++]
        l = line.split(" ")
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
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            if (j > 0)
                stream.write(" ")
            stream.write(`${a[i][j]}`)
        }
        stream.write("\n")
    }
}

a = matrix_read("/dev/stdin")
M = a.length
N = a[0].length
matrix_write(process.stdout, M, N, a)
