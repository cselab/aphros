include(mh.m4)dnl
#!/usr/bin/env node
"use strict"
mh_include(matrix.js)dnl
mh_include(partstr.js)dnl

var argv, nh, hl, hp, eta, a, t, file, u, M, N
var ends, p, end, partstr, i, j, ne, n

argv = process.argv

nh = 2
hl = 2
hp = 4 / (2.0*nh)
eta = 0.5
a = 0.0
t = 0.0
n = 2*nh + 1

argv.shift(); argv.shift()
file = argv.shift()
u = matrix_read(file)
M = u.length
N = u[0].length
matrix_halo_zero(M, N, hl, u)
ends = matrix_new(M, N)
partstr_vof_ends(M, N, u, /**/ ends)

partstr = new Partstr(nh, hl, hp, eta)

i = 1; j = 2
end = partstr_cell_ends(M, N, i, j, ends)
ne = end.length
p = [i + 0.5, j + 0.5]
partstr.start(ne, end, a, t, p)
for (i = 0; i < 10000; i++)
    partstr.step()

partstr_cell_ends_gnuplot(process.stdout, ne, partstr.ends)
process.stdout.write("\n\n")
partstr_force_write(process.stdout, n, partstr.xx, partstr.ff)
