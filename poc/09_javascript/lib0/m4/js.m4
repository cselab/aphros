divert(-1)dnl
include(mh.m4)dnl
mh_define(<[[FU]]>, <[[function $1($2) {
    if (arguments.length != $1.length)
       throw new Error(`$1: wrong numer of arguments\n`)dnl
$3dnl
}]]>)dnl
mh_define(<[[INT]]>, <[[dnl
if (!Number.isInteger($1)) throw new Error(`$1 = ${$1}\n`)]]>)dnl
divert<[[]]>dnl
