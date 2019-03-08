#!/usr/bin/env gnuplot

set macros
reset

set terminal pdfcairo color font "Helvetica,14" size 3.2,2.2 enhanced

c0='#1f77b4'
c1='#ff7f0e'
c2='#2ca02c'
c3='#d62728'
c4='#9467bd'

set linetype  1 lc rgb c0 lw 2
set linetype  2 lc rgb c2 lw 2
set linetype  3 lc rgb c1 lw 2

set style line 1 lt 1 pt 7 ps 0.5
set style line 2 lt 2 pt 5 ps 0.5
set style line 3 lt 3 pt 8 ps 0.5

set key bottom

set output "a.pdf"

# frame index to simulation time
sim = 0.01
# bubble radius
r = 0.15
# simulation time to dimensionless time
tau = 2 ** 1.5

set xrange [0:0.4]

plot "ref/rn" w l ls 1 t 'exp' , \
  "neck" u ($0*sim*tau):($2/r-1+0.03) w lp ls 2 t '512' , \
  "neck256" u ($0*sim*tau):($2/r-1+0.043) w lp ls 3 t '256' , \

