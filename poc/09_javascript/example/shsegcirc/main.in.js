undivert(matrix.js)dnl
undivert(partstr.js)dnl
changequote()

const stdout = process.stdout
const stderr = process.stdout

var j = 2, k, l, d
var k = parseFloat(process.argv[j++])

var ax = parseFloat(process.argv[j++])
var ay = parseFloat(process.argv[j++])

var bx = parseFloat(process.argv[j++])
var by = parseFloat(process.argv[j++])

var xx = parseFloat(process.argv[j++])
var xy = parseFloat(process.argv[j++])

stdout.write(`${partstr_shsegcirc(k, [ax, ay], [bx, by], [xx, xy]o)}\n`)
