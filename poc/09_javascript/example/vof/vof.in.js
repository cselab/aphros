include(mh.m4)dnl
#!/usr/bin/env node
"use strict"
mh_include(vof.js)dnl

function f(x, y) {
    return x > 0
}

var vof = new Vof(0.1, 0.1, f)

console.log(vof.cell(1, 1))
