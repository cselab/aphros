undivert(partstr.js)dnl
changequote()

const AX = 0, AY = 1, BX = 2, BY = 3
var nx, ny, a
var seg = Array(4)

partstr_ends(nx, ny, a, seg)
process.stdout.write(`${seg[AX]} ${seg[AY]} ${seg[BX]} ${seg[BY]}\n`)
