include(mh.m4)dnl
#!/usr/bin/env node
"use strict"
mh_include(matrix.js)dnl
mh_include(vof.js)dnl

function f(x, y, param) {
    var x0, y0, r
    x0 = param.x0
    y0 = param.y0
    r = param.r
    x -= x0
    y -= y0
    return x*x + y*y < r*r
}

var u, lx, ly, dx, dy, M, N, i, j, x0, y0, r
var vof, param

M = N = 5

lx = ly = 1

dx = lx/M
dy = ly/N

x0 = y0 = 0.5; r = 0.33
param = {x0: x0, y0: y0, r: r}

vof = new Vof(dx, dy, f, param)
u = matrix_new(M, N)
vof.grid(M, N, u)

matrix_write(process.stdout, M, N, u)
