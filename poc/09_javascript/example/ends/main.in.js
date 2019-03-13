undivert(partstr.js)dnl
changequote()

const stdout = process.stdout
const stderr = process.stderr

const AX = 0, AY = 1, BX = 2, BY = 3
var nx, ny, a
var seg = Array(4)

i = 2
nx = parseFloat(process.argv[i++])
ny = parseFloat(process.argv[i++])
a  = parseFloat(process.argv[i++])

partstr_ends(nx, ny, a, seg)
stdout.write(`${seg[AX]} ${seg[AY]} ${seg[BX]} ${seg[BY]}\n`)
