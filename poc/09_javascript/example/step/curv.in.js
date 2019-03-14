include(mh.m4)dnl
#!/usr/bin/node
mh_include(matrix.js)dnl
mh_include(partstr.js)dnl

const stdout = process.stdout
const stderr = process.stderr
var argv = process.argv

argv.shift(); argv.shift()
var t = parseFloat(argv.shift())
var hp = 2.0

var k = partstr_curv(hp, t)

stdout.write(`${k}\n`)
