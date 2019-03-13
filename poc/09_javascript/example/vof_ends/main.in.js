undivert(matrix.js)dnl
undivert(partstr.js)dnl
changequote()

const AX = 0, AY = 1, BX = 2, BY = 3
const stdout = process.stdout
const stderr = process.stderr
const hl = 2

var i = 2
var file = process.argv[i++]
var u = matrix_read(file)
var M = u.length
var N = u[0].length

matrix_halo_zero(M, N, hl, u)

var ends = matrix_new(M, N)
partstr_vof_ends(M, N, u, ends)

var i, j, eb
for (i = 0; i < M; i++)
    for (j = 0; j < M; j++) {
        e = ends[i][j]
        stdout.write(`${e[AX]} ${e[AY]}\n`)
        stdout.write(`${e[BX]} ${e[BY]}\n`)
        stdout.write("\n")
    }
