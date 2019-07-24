#!/usr/bin/env gnuplot

set macros
reset

set terminal pdfcairo color font "Helvetica,14" size 3.2,2.2 enhanced


c0='#1f77b4'
c1='#ff7f0e'
c2='#2ca02c'
c3='#d62728'
c4='#9467bd'

set yrange [1e-16:1]

set linetype  1 lc rgb c0 lw 3 
set linetype  2 lc rgb c2 lw 3
set linetype  3 lc rgb c1 lw 3
set linetype  4 lc rgb c3 lw 3
set linetype  5 lc black lw 3 dt 2
set linetype cycle 4

set style line 1 lt 1 pt 7 ps 0.5 
set style line 2 lt 2 pt 5 ps 0.5 
set style line 3 lt 3 pt 8 ps 0.5 

set logscale y

set format y "10^{%L}"
set output "a.pdf"

set xtics 20

set xlabel "m"
set ylabel "E_m"

#set nokey

#Shadecolor = "#80E0A080"

#set arrow from 20, graph 0 to 20, graph 1 nohead lt 5

ll = system("echo iter*")

plot for [f in ll] f u 0:1 every ::1 w l t f , \
  for [f in ll] f u 0:2:3 every ::1 with filledcurve fs transparent solid 0.5 t '' lt 1 ,  \
  1e-5 w l lt 5

