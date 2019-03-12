#!/usr/bin/env node
undivert(partstr.js)dnl
changequote()dnl
var stdout = process.stdout
var stderr = process.stdout

const AX = 0, AY = 1, BX = 2, BY = 3
var nx, ny, a
var seg = Array(4)

var i = 2
var nx = parseFloat(process.argv[i++])
var ny = parseFloat(process.argv[i++])
var u  = parseFloat(process.argv[i++])

var a = partstr_line(nx, ny, u)
stdout.write(`${a}\n`)
