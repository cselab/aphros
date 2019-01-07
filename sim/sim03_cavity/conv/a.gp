set macros
reset

set terminal pdfcairo color font "Helvetica,14" size 4,2.2 enhanced

c0='#1f77b4'
c1='#ff7f0e'
c2='#2ca02c'
c3='#d62728'
c4='#9467bd'

set linetype  1 lc rgb c0 lw 2
set linetype  2 lc rgb c2 lw 2
set linetype  3 lc rgb c1 lw 2
set linetype  4 lc rgb c3 lw 2
set linetype  5 lc rgb "black" lw 1 dt 1
set linetype cycle 4

set style line 1 lt 1 pt 7 ps 0.
set style line 2 lt 2 pt 5 ps 0. pi 2
set style line 3 lt 3 pt 7 ps 0. pi 3
set style line 4 lt 4 pt 2 ps 0. pi 3

set key vertical outside spacing 1.25 

first(x) = ($0 > 0 ? base : base = x)

ss="ch064 ch128 ch256 ch512"
tt="64 128 256 512"

set output "vy.pdf"
set xtics 0.2
plot for [i=1:words(ss)] "<paste ".word(ss,i)."/{x,vy}" w lp ls i t word(tt,i)

set output "dvy.pdf"
set logscale y 10
set xtics auto
set ytics auto
plot for [i=1:words(ss)] word(ss,i)."/dvy" u 1:(abs($2)) w lp ls i t word(tt,i)

set output "er.pdf"
set logscale xy 2
plot \
  "er" u "nx":(first(column("e2")), column("e2")/base) w lp ps 0.5 pt 7 t "L2" \
  , "er" u "nx":(first(column("em")), column("em")/base) w lp ps 0.5 pt 7 t "max" \
  , (64/x) dt 4 t 'first' \
  , (64/x)**(2) dt 2 t 'second' \



