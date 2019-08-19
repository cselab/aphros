#!/usr/bin/env gnuplot

set macros
#reset

set terminal pdfcairo color font "Helvetica,14" size 4.2,2.2 enhanced

#set key reverse bottom right Left
set key outside

c0='#1f77b4'
c1='#ff7f0e'
c2='#2ca02c'
c3='#d62728'
c4='#9467bd'
c5='#8c564b'

if (!exists("ss")) {
  ss = "bap ba"
}

ss = " \
univel_dim2_dom1 \
univel_dim2_dom100 \
univel_dim2_dom1_nofix \
univel_dim2_dom2_nofix \
"

tt = "\
  fix,L=1 \
  fix,L=100 \
  orig,L=1 \
  orig,L=2 \
  "

set linetype  1 lc rgb c0 lw 2
set linetype  2 lc rgb c2 lw 2
set linetype  3 lc rgb c1 lw 2
set linetype  4 lc rgb c3 lw 2
set linetype  5 lc rgb c4 lw 2
set linetype  6 lc rgb c5 lw 2
set linetype  9 lc rgb "black" lw 2


set linetype cycle 4

set style line 1 lt 1 pt 7 ps 0.5 pi 1
set style line 2 lt 1 pt 5 dt '--' ps 0.5 pi 1
set style line 3 lt 3 pt 8 ps 0.5 pi 1
set style line 4 lt 3 pt 2 dt '--' ps 0.5 pi 1
set style line 5 lt 5 pt 2 ps 0.5 pi 1
set style line 6 lt 6 pt 2 ps 0.5 pi 1


p="/erba"
m='plot for [i=1:words(ss)] word(ss,i).p u "cpr":v w l ls i t word(tt,i)'
s='set output "a".v.".pdf"' 
set logscale x 2
set xrange [0.5:32]
set xlabel "R / h"

set format y "10^{%L}"
set yrange [1e-8:1]
set logscale y
#set ylabel "We_{rms}"
#v="vl2x" ; @s ; @m
set ylabel "We_{max}" offset -1
v="vlmx" ; @s ; @m 

unset format y
unset logscale y
set yrange [0:2]
set ylabel "{/Symbol D}p / p_L"
v="pd" ; @s ; @m

