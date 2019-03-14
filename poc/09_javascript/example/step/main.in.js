include(mh.m4)dnl
#!/usr/bin/node
mh_include(matrix.js)dnl
mh_include(partstr.js)dnl

const X = 0, Y = 1
const stdout = process.stdout
const stderr = process.stderr
const argv = process.argv
var ename, pname, k, ends, xx, ne, np, eta, ff, a, t

argv.shift(); argv.shift()
ename = argv.shift()
k = parseFloat(argv.shift())
a = 0.1
t = 0.2

nh = 2
np = 2*nh + 1
hp = 1.0
p = [0, 0]

ends = partstr_ends_read(ename)
ne = ends.length

xx = matrix_new(2*nh + 1, 2)
partstr_part(nh, hp, p, a, t, /**/ xx)

eta = 1
ff = matrix_new(np, 2)
k = partstr_curv(hp, t)
_msg(`ff: ${ne} ${np} ${k} ${eta}\n`)
partstr_force(ne, ends, np, xx, k, eta,     ff)

_msg(`ff: ${ff}\n`)

var p = [0, 0]
var State = { p: p, a: a, t: t }
var hp = 1
partstr_step(nh, ff, eta, hp, State)
