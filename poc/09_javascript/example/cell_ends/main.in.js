undivert(matrix.js)dnl
undivert(partstr.js)dnl
changequote()

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
