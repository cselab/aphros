include(mh.m4)dnl
#!/usr/bin/env node
mh_include(matrix.js)dnl
mh_include(partstr.js)dnl

function Partstr(nh, hl, hp, eta) {
    "use strict"
    this.nh = nh
    this.hl = hl
    this.hp = hp
    this.eta = eta
    this.new = function(M, N, ends, i, j, a, t) {
        if (!Array.isArray(ends))
            throw new Error(`ends is not an array: ${ends}`)
        this.ends = partstr_cell_ends(M, N, i, j, ends)
        this.a = a
        this.t = t
        this.p = [i + 0.5, j + 0.5]
    }
    this.step = function() {
        var nh, p, a, t, hp, eta, ends, np, ff, State, k, ne, xx
        nh = this.nh
        p = this.p
        a = this.a
        t = this.t
        hp = this.hp
        eta = this.eta
        ends = this.ends
        
        np = 2*nh + 1
        ff = matrix_new(np, 2)
        State = { p: p, a: a, t: t }
        k = partstr_curv(hp, t)
        ne = ends.length
        xx = matrix_new(np, 2)
        partstr_part(nh, hp, p, a, t, /**/ xx)
        partstr_force(ne, ends, np, xx, k, eta,     ff)
        partstr_step(nh, ff, this.eta, this.hp, /*io*/ State)
        this.a = State.a
        this.t = State.t
        this.p = State.p
        this.k = k
        this.xx = xx
        this.ff = ff
    }
}

argv = process.argv

nh = 2
hl = 2
hp = 4 / (2.0*nh)
eta = 1.0
a = 0.0
t = 0.0

argv.shift(); argv.shift()
file = argv.shift()
u = matrix_read(file)
M = u.length
N = u[0].length
matrix_halo_zero(M, N, hl, u)
ends = matrix_new(M, N)
partstr_vof_ends(M, N, u, /**/ ends)

partstr = new Partstr(nh, hl, hp, eta)
i = 1; j = 1

partstr.new(M, N, ends, i, j, a, t)
for (i = 0; i < 2; i++)
    partstr.step()

ne = partstr.ends.length
partstr_cell_ends_gnuplot_write(process.stdout, ne, partstr.ends)
process.stdout.write("\n\n")
partstr_force_write(process.stdout, 2*nh + 1, partstr.xx, partstr.ff)
