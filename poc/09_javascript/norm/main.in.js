undivert(lib.js)dnl
undivert(matrix.js)dnl
undivert(partstr.js)dnl
changequote()

var stdout = process.stdout
var n = 5
var h = 2
var u = matrix_new(n, n)
InitGrid(n, n, u)
matrix_halo_zero(n, n, h, u)
matrix_lh_write(stdout, -h, -h, n + h, n + h, u)

var p = partstr_norm(2, 0, u)
stdout.write(`norm: ${p}\n`)

const X = 0, Y = 1
a = partstr_line(p[X], p[Y], u[2][0])
stdout.write(`a: ${a}\n`)
