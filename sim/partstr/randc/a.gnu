#!/usr/bin/env gnuplot

set terminal pdfcairo color font "Helvetica,14" size 3.2,2.2 enhanced
set macro
reset

set nokey

c0='#1f77b4'
c1='#ff7f0e'
c2='#2ca02c'
c3='#d62728'
c4='#9467bd'

set linetype  1 lc rgb c0 lw 2
set linetype  2 lc rgb c2 lw 2
set linetype  3 lc rgb c1 lw 2
set linetype  4 lc rgb c3 lw 2
set linetype  5 lc rgb c4 lw 2
set linetype  6 lc rgb "black" lw 2
set linetype cycle 6

set style line 1 lt 1 pt 7 ps 0.5
set style line 2 lt 2 pt 5 ps 0.5
set style line 3 lt 3 pt 8 ps 0.5
set style line 4 lt 4 pt 2 ps 0.5

ss="ch ba bap"

xx='(column("cpr"))'

s='pt 7 ps 0.5'
m='for [i=1:words(ss)] word(ss,i)."/kerravg" \
    u @xx:@v w lp t word(ss,i) ls i'
mfs='for [i=1:words(ss)] word(ss,i)."/kerravg" \
    u @xx:@v."sl":@v."sh" w filledcurve t "" ls i fs transparent solid 0.5'
mf='for [i=1:words(ss)] word(ss,i)."/kerravg" \
    u @xx:@v."l":@v."h" w filledcurve t "" ls i fs transparent solid 0.25'

set logscale x 2

set xrange [0.5:32]
set xlabel "R / h"

set output "aem.pdf"
set logscale y
set ylabel "L_{/Symbol \245}"
set yrange [0.0001:10]
v='"em"'
plot \
1/x**2 dt '-.' lc "black" lt 5 , \
@m , @mfs \
#, "mdhf/lm" u ($1/2):"e" t "mdhf" w lp ls 3 \
#, "ref/ba" u ($1/2):"max" t "baweb" w lp ls 4 \

set output "ae2.pdf"
set logscale y
set ylabel "L_2"
set yrange [0.0001:10]
v='"em"'
v='"e2"'
plot \
1/x**2 dt '-.' lc "black" lt 5 , \
@m , @mfs  \
#, "mdhf/l2" u ($1/2):"e" t "mdhf" w lp ls 3 \
#, "ref/ba" u ($1/2):"rms" t "baweb" w lp ls 4 \

set output "avf.pdf"
unset logscale y
set yrange [0:3]
set ylabel "V / V_{exact}"
v='"vf"'
plot @m , @mfs
