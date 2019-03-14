include(mh.m4)dnl
#!/usr/bin/node
mh_include(matrix.js)dnl
mh_include(partstr.js)dnl

const X = 0, Y = 1
const stdout = process.stdout
const stderr = process.stderr
const argv = process.argv
var ename, pname, k, ends, xx, ne, np, eta, ff

argv.shift(); argv.shift()

ename = argv.shift()
pname = argv.shift()
k = parseFloat(argv.shift())
ends = partstr_ends_read(ename)
xx   = matrix_read(pname)

ne = ends.length
np = xx.length

eta = 1
ff = matrix_new(np, 2)
partstr_force(ne, ends, np, xx, k, eta, /**/ ff)

stderr.write(`ne = ${ne}\n`)

//partstr_ends_gnuplot_write(stdout, ne, ends)
partstr_force_write(stdout, np, xx, ff)
