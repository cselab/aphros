undivert(matrix.js)dnl
undivert(partstr.js)dnl
undefine(shift)dnl
changequote()

const stdout = process.stdout
const stderr = process.stdout
const argv = process.argv

argv.shift(); argv.shift()
var ax = parseFloat(argv.shift())
var ay = parseFloat(argv.shift())

var bx = parseFloat(argv.shift())
var by = parseFloat(argv.shift())

var xx = parseFloat(argv.shift())
var xy = parseFloat(argv.shift())

stdout.write(`${ partstr_nearest([ax, ay], [bx, by], [xx, xy])}\n`)
