include(mh.m4)dnl
#!/usr/bin/node
mh_include(matrix.js)dnl
mh_include(partstr.js)dnl

const X = 0, Y = 1
const stdout = process.stdout
const stderr = process.stderr
const argv = process.argv
var ename, pname, k, ends, xx, ne, np, eta, ff, a, t, State

argv.shift(); argv.shift()
ename = argv.shift()
a = 0.2
t = 0.3
p = [0.1, 0.2]
nh = 2
np = 2*nh + 1
hp = 1.0
eta = 2.0

ends = partstr_ends_read(ename)
ne = ends.length

xx = matrix_new(np, 2)
partstr_part(nh, hp, p, a, t, /**/ xx)

ff = matrix_new(np, 2)
k = partstr_curv(hp, t)
partstr_force(ne, ends, np, xx, k, eta,     ff)

State = { p: p, a: a, t: t }
partstr_step(nh, ff, eta, hp, /*io*/ State)
_msg(`${State.p} ${State.a} ${State.t}\n`)
