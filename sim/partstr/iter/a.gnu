#!/usr/bin/env gnuplot

set macros
reset

set terminal pdfcairo color font "Helvetica,14" size 3.2,2.2 enhanced


c0='#1f77b4'
c1='#ff7f0e'
c2='#2ca02c'
c3='#d62728'
c4='#9467bd'

set linetype  1 lc rgb c0 lw 3 
set linetype  2 lc rgb c2 lw 3
set linetype  3 lc rgb c1 lw 3
set linetype  4 lc rgb c3 lw 3
set linetype cycle 4

set style line 1 lt 1 pt 7 ps 0.5 
set style line 2 lt 2 pt 5 ps 0.5 
set style line 3 lt 3 pt 8 ps 0.5 

set logscale y

set output "a.pdf"

Shadecolor = "#80E0A080"
#
plot 'iter.dat' u 0:1 w l t '', \
  '' using 0:2:3 with filledcurve fc rgb Shadecolor t ''

