include(mh.m4)dnl
#!/usr/bin/node
mh_include(matrix.js)dnl
mh_include(partstr.js)dnl

const stderr = process.stderr
const stdout = process.stdout
function msg(s) { return stderr.write(s) }

const nh = 2, hp = 1, a = 0.1, t = 0.3

var n = 2*nh + 1

var da = matrix_new(n, 2)
var dt = matrix_new(n, 2)

partstr_dxda(nh, hp, a, t, /**/ da)
partstr_dxdt(nh, hp, a, t, /**/ dt)

matrix_write(stdout, n, 2, da)
matrix_write(stdout, n, 2, dt)
