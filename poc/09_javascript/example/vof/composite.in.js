include(mh.m4)dnl
#!/usr/bin/env node
"use strict"
mh_include(matrix.js)dnl
mh_include(vof.js)dnl

var u, w, lx, ly, dx, dy, M, N, i, j, x0, y0, a, b
var vof, Param, f

M = N = 5

lx = ly = 1

dx = lx/M
dy = ly/N

a = 0.35; b = 0.20
Param = {}
Param.param = [{x0: 0.4, y0: 0.5, a: a, b: b},
               {x0: 0.3, y0: 0.5, a: a, b: b}]
Param.f = [vof_ellipse, vof_ellipse]

vof = new Vof(dx, dy, vof_comosite, Param)
u = matrix_new(M, N)
vof.grid(M, N, u)

w = matrix_transpose(M, N, u)
matrix_write(process.stdout, N, M, w)
