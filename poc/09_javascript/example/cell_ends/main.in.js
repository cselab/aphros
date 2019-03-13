undivert(matrix.js)dnl
undivert(partstr.js)dnl
changequote()

const AX = 0, AY = 1, BX = 2, BY = 3
const stdout = process.stdout
const stderr = process.stderr

const hl = 2

var k = 2, i, j, e
var file = process.argv[k++]
var i = parseInt(process.argv[k++])
var j = parseInt(process.argv[k++])
var u = matrix_read(file)
var M = u.length
var N = u[0].length
matrix_halo_zero(M, N, hl, u)
var ends = matrix_new(M, N)

partstr_vof_ends(M, N, u, ends)
e = partstr_cell_ends(M, N, i, j, ends)
partstr_ends_write(stdout, e)
