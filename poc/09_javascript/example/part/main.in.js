include(mh.m4)dnl
undivert(matrix.js)dnl
undivert(partstr.js)dnl
changequote()

const stdout = process.stdout
const stderr = process.stderr

var k = 2, i, j, e
var nh = parseInt(process.argv[k++])
var hp = parseFloat(process.argv[k++])
var a = parseFloat(process.argv[k++])
var t = parseFloat(process.argv[k++])

var p = [0, 0]

var xx = partstr_part(nh, hp, p, a, t)
var M = xx.length
matrix_write(stdout, M, 2, xx)
