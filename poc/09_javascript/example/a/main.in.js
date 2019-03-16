include(mh.m4)dnl
#!/usr/bin/env node
"use strict"
mh_include(matrix.js)dnl
mh_include(partstr.js)dnl

var argv, nh, hl, hp, eta, eps, a, t, file, u, M, N, e
var ends, p, end, partstr, i, itermax, j, ne, n, eps, k
var AX = 0, AY = 1, BX = 2, BY = 3

argv = process.argv

nh = 4
hp = 4 / (2.0*nh)
eta = 0.5
t = 0.0
n = 2*nh + 1

argv.shift(); argv.shift()
file = argv.shift()
end = partstr_ends_read(file)
partstr = new Partstr(nh, hp, eta)

e = end[0]
p = [0, 1]
a = Math.PI

ne = end.length
partstr.start(ne, end, a, t, p)
eps = 1e-5
itermax = 20
partstr.converge(1e-5, 20)

partstr_cell_ends_gnuplot(process.stdout, ne, partstr.ends)
process.stdout.write("\n\n")
partstr_force_write(process.stdout, n, partstr.xx, partstr.ff)
