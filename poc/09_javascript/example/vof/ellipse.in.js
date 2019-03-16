include(mh.m4)dnl
#!/usr/bin/env node
"use strict"
mh_include(matrix.js)dnl
mh_include(vof.js)dnl

function f(x, y, param) {
    var x0, y0, a, b
    x0 = param.x0
    y0 = param.y0

    a = param.a
    b = param.b

    x -= x0
    y -= y0
    return (x*x)/(a*a) + (y*y)/(b*b) < 1
}

var u, w, lx, ly, dx, dy, M, N, i, j, x0, y0, a, b
var vof, param

M = N = 5

lx = ly = 1

dx = lx/M
dy = ly/N

x0 = y0 = 0.5; a = 0.35; b = 0.20
param = {x0: x0, y0: y0, a: a, b: b}

vof = new Vof(dx, dy, f, param)
u = matrix_new(M, N)
vof.grid(M, N, u)

w = matrix_transpose(M, N, u)
matrix_write(process.stdout, N, M, w)
