include(mh.m4)dnl
#!/usr/bin/env node
"use strict"
mh_include(vof.js)dnl

function f(x, y) { return x > 0 && y > 0; }
console.log(vof_trapz(100, 100, f, -0.5, -0.5, 0.5, 0.5))

