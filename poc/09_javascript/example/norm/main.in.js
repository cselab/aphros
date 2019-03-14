include(mh.m4)dnl
mh_include(lib.js)dnl
mh_include(matrix.js)dnl
mh_include(partstr.js)dnl

var stdout = process.stdout
var n = 5
var h = 2
var u = matrix_new(n, n)
var a = matrix_new(n, n)

InitGrid(n, n, u)
matrix_halo_zero(n, n, h, u)
partstr_vof_line(n, n, u, a)
matrix_write(stdout, n, n, a)
