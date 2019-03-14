include(mh.m4)dnl
mh_include(matrix.js)dnl
mh_include(partstr.js)dnl

const stdout = process.stdout
const stderr = process.stderr

var j = 2, k, l, d
var k = parseFloat(process.argv[j++])

var ax = parseFloat(process.argv[j++])
var ay = parseFloat(process.argv[j++])

var bx = parseFloat(process.argv[j++])
var by = parseFloat(process.argv[j++])

var xx = parseFloat(process.argv[j++])
var xy = parseFloat(process.argv[j++])

stdout.write(`${partstr_shsegcirc(k, [ax, ay], [bx, by], [xx, xy]o)}\n`)
