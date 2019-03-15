include(mh.m4)dnl
#!/usr/bin/env node
"use strict"
mh_include(matrix.js)dnl
mh_include(partstr.js)dnl

function Partstr(nh, hl, hp, eta, M, N) {
    const n = 2*nh + 1
    this.nh = nh
    this.hl = hl
    this.hp = hp
    this.eta = eta
    this.ff = matrix_new(n, 2)
    this.xx = matrix_new(n, 2)
    this.M = M
    this.N = N
    this.new = function(ends, i, j, a, t) {
        if (!Array.isArray(ends))
            throw new Error(`ends is not an array: ${ends}`)
        this.ends = partstr_cell_ends(this.M, this.N, i, j, ends)
        this.a = a
        this.t = t
        this.p = [i + 0.5, j + 0.5]
    }

    this.step = function() {
        var nh, p, a, t, hp, eta, ends, ff, State, k, ne, xx
        nh = this.nh
        const n = 2*nh + 1
        hp = this.hp
        eta = this.eta
        ends = this.ends
        xx = this.xx
        ff = this.ff
        p = this.p; a = this.a; t = this.t
        k = partstr_curv(hp, t)
        ne = ends.length
        partstr_part(nh, hp, p, a, t, /**/ xx)
        partstr_force(ne, ends, n, xx, k, eta,     ff)
        State = { p: p, a: a, t: t }
        partstr_step(nh, ff, eta, hp, /*io*/ State)
        this.a = State.a; this.t = State.t; this.p = State.p
        this.k = k
        this.xx = xx
        this.ff = ff
    }
}

var argv, nh, hl, hp, eta, a, t, file, u, M, N, ends, partstr, i, j, ne, n

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

partstr = new Partstr(nh, hl, hp, eta, M, N)
i = 2; j = 4
partstr.new(ends, i, j, a, t)
for (i = 0; i < 10000; i++)
    partstr.step()

ne = partstr.ends.length
partstr_cell_ends_gnuplot_write(process.stdout, ne, partstr.ends)
process.stdout.write("\n\n")
partstr_force_write(process.stdout, n, partstr.xx, partstr.ff)
