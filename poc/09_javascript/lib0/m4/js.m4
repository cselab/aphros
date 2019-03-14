include(mh.m4)dnl
mh_define(<[[FU]]>, <[[function $1($2) {
    if (arguments.length != $1.length)
       throw new Error(`$1: wrong numer of arguments\n`)dnl
$3dnl
}]]>)dnl
dnl
FU(f, <[[x, y]]>, <[[
    console.log(`${x}`)
    console.log(`${y}`)
]]>)
