include(mh.m4)dnl
#!/usr/bin/env node
"use strict"
mh_include(matrix.js)dnl
mh_include(vof.js)dnl

var u, w, lx, ly, dx, dy, M, N, i, j, x0, y0, a, b
var vof, param, f

M = N = 5

lx = ly = 1

dx = lx/M
dy = ly/N

x0 = y0 = 0.5; a = 0.35; b = 0.20
param = {x0: x0, y0: y0, a: a, b: b}
f = vof_ellipse

vof = new Vof(dx, dy, f, param)
u = matrix_new(M, N)
vof.grid(M, N, u)

w = matrix_transpose(M, N, u)
matrix_write(process.stdout, N, M, w)
