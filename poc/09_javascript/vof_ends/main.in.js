undivert(matrix.js)dnl
undivert(partstr.js)dnl
changequote()

var stdout = process.stdout
var stderr = process.stdout

var i = 2
var file = process.argv[i++]
var u = matrix_read(file)

var M = u.length
var N = u[0].length

var hl = 2
matrix_halo_zero(M, N, hl, u)

var ends = matrix_new(M, N)
partstr_vof_ends(M, N, u, ends)

matrix_write(stdout, M, N, ends);
