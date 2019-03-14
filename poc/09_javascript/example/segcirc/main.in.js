include(mh.m4)dnl
mh_include(matrix.js)dnl
mh_include(partstr.js)dnl

const stdout = process.stdout
const stderr = process.stderr

var j = 2, k, l, d
var k = parseFloat(process.argv[j++])
var l = parseFloat(process.argv[j++])
var d = parseFloat(process.argv[j++])

stdout.write(`${partstr_segcirc(k, l, d)}\n`)
