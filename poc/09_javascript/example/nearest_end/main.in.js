include(mh.m4)dnl
undivert(matrix.js)dnl
undivert(partstr.js)dnl
undefine(shift)dnl
changequote()

const X = 0, Y = 1
const stdout = process.stdout
const stderr = process.stderr
const argv = process.argv
var ax, ay, bx, by, xx, xy, y
var k, ends, n, file

argv.shift(); argv.shift()
file = argv.shift()
k = parseFloat(argv.shift())
xx = parseFloat(argv.shift())
xy = parseFloat(argv.shift())

ends = partstr_ends_read(file)
n = ends.length
y = partstr_nearest_ends(n, ends, [xx, xy], k)

stderr.write(`${y[X]} ${y[Y]}\n`)
partstr_ends_gnuplot_write(stdout, ends)
