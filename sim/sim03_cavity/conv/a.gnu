#!/usr/bin/env gnuplot

set macros
reset

set terminal pdfcairo color font "Helvetica,14" size 4,2.2 enhanced

c0='#1f77b4'
c1='#ff7f0e'
c2='#2ca02c'
c3='#d62728'
c4='#9467bd'

set linetype  1 lc rgb c0 lw 1
set linetype  2 lc rgb c2 lw 1
set linetype  3 lc rgb c1 lw 1
set linetype  4 lc rgb c3 lw 1
set linetype  5 lc rgb c4 lw 1
set linetype  9 lc rgb "black" lw 1 dt 1
set linetype cycle 5

set style line 1 lt 1 pt 7 ps 0.
set style line 2 lt 2 pt 5 ps 0. pi 2
set style line 3 lt 3 pt 7 ps 0. pi 3
set style line 4 lt 4 pt 2 ps 0. pi 3
set style line 5 lt 5 pt 2 ps 0. pi 3
set style line 9 lt 5 pt 0 ps 0. pi 3

set key vertical outside spacing 1.25 

first(x) = ($0 > 0 ? base : base = x)

ss="ch064 ch128 ch256 ch512 ch1024"
tt="64 128 256 512 1024"

set output "vy.pdf"
set xtics 0.2
plot for [i=1:words(ss)] "<paste ".word(ss,i)."/{x,vy}" w lp ls i t word(tt,i)

set output "dvy.pdf"
set logscale y 10
set xtics auto
set ytics auto
plot for [i=1:words(ss)-1] word(ss,i)."/dvy" u 1:(abs($2)) w lp ls i t word(tt,i)

f1(x) = a1 / x
f2(x) = a2 / x ** 2
set fit logfile "/dev/null"
fit log(f1(x)) 'er' using "nx":(log(column("e2"))) via a1
fit log(f2(x)) 'er' using "nx":(log(column("e2"))) via a2

set output "er.pdf"
set logscale xy 2
plot \
  "er" u "nx":"e2" w lp ps 0.5 pt 7 t "L2" \
  , "er" u "nx":"em" w lp ps 0.5 pt 7 t "max" \
  , f1(x) dt 4 t 'first' \
  , f2(x) dt 2 t 'second' \

