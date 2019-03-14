include(mh.m4)dnl
#!/usr/bin/node
mh_include(matrix.js)dnl
mh_include(partstr.js)dnl

const stdout = process.stdout

const nh = 2, a = 0.1, t = 0.3
var xx = matrix_new(2*nh + 1, 2)

partstr_dxdt(nh, a, t, /**/ xx)

stdout.write("Hello world!\n")
