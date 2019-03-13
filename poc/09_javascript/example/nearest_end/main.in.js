undivert(matrix.js)dnl
undivert(partstr.js)dnl
undefine(shift)dnl
changequote()

const stdout = process.stdout
const stderr = process.stdout
const argv = process.argv

argv.shift(); argv.shift()
var ax, ay, bx, by, xx, xy
ax = parseFloat(argv.shift())
ay = parseFloat(argv.shift())

bx = parseFloat(argv.shift())
by = parseFloat(argv.shift())

xx = parseFloat(argv.shift())
xy = parseFloat(argv.shift())

stdout.write(`${ partstr_nearest([ax, ay], [bx, by], [xx, xy])}\n`)
