#!/usr/bin/env gnuplot

#!rm -f {vx,vy}_*.{xmf,h5}
#!./rs

termsize="3,2"
load "style.gp"

set key font ",8"

flt(x) = system("echo '".x."' | tr '_' '-'")

set key bottom

set xlabel lbt
set xtics 5
set xrange [0:40]

nn="-E '51'"

llstat=system("echo ref2/quick_ext40_lim01/stat_{re40_nd026,re40_nd051}")
llsep=system("echo ref2/quick_ext40_lim01/sep_{re40_nd026,re40_nd051}")

set output "drag.pdf"
set ylabel "{/:Italic C}_D"
plot [][1.2:2.5] for [f in llstat]\
  f u "t":(column("dragx")*2) w l t flt(f) ,\
1.52 w l ls 1 dt '-' t "Choi2007,Re=40" ,\
2.02 w l ls 2 dt '-' t "Choi2007,Re=20"

set output "vdrag.pdf"
set ylabel "{/:Italic C}_{D,vis}"
plot [][0.45:0.8] for [f in llsep]\
  f u "t":(column("vdragx")*2) w l t flt(f) ,\

set output "sep.pdf"
set ylabel "{/:Italic L} / {/:Italic D}"
plot [][0:2.5] for [f in system(\
  "echo ref*/quick*/sep_* | tr ' ' '\n'  | grep lim01 | grep ".nn)." ".system("echo sep_*_*")]\
  f w l t flt(f) ,\
2.25 w l lw 2 lc rgb "black" t "Choi2007,Re=40" ,\
0.90 w l lw 2 lc rgb "black" t "Choi2007,Re=20"
