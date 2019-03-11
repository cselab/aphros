undivert(lib.js)dnl
undivert(matrix.js)dnl
changequote()

var stdout = process.stdout
var u = []
InitGrid(5, 5, u)
matrix_lh_write(stdout, -3, -3, 5, 5, u)
