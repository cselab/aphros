undivert(lib.js)dnl
undivert(matrix.js)dnl
changequote()

var u = []
InitGrid(5, 5, u)
matrix_write(process.stdout, 5, 5, u)
