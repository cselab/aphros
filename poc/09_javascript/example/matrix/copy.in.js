include(mh.m4)dnl
#!/usr/bin/env node
"use strict"
mh_include(matrix.js)dnl

var a, b, M, N

a = matrix_read(process.argv[2])
M = a.length
N = a[0].length

b = matrix_copy(M, N, a)

matrix_write(process.stdout, M, N, b)
