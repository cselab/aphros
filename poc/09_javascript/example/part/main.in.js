include(mh.m4)dnl
#!/usr/bin/env node
mh_include(matrix.js)dnl
mh_include(partstr.js)dnl

const stdout = process.stdout
const stderr = process.stderr

var k = 2, i, j, e
var nh = parseInt(process.argv[k++])
var hp = parseFloat(process.argv[k++])
var a = parseFloat(process.argv[k++])
var t = parseFloat(process.argv[k++])

var p = [0, 0]

var xx = matrix_new(2*nh + 1, 2)
partstr_part(nh, hp, p, a, t, /**/ xx)
var M = xx.length
matrix_write(stdout, M, 2, xx)
