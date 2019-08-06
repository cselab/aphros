#!/usr/bin/env gnuplot

set macros
reset

set terminal pdfcairo color font "Helvetica,14" size 3.2,2.2 enhanced

#set nokey
c0='#1f77b4'
c1='#ff7f0e'
c2='#2ca02c'
c3='#d62728'
c4='#9467bd'

if (!exists("ss")) {
  ss = "bap ba"
}

set linetype  1 lc rgb c0 lw 2
set linetype  2 lc rgb c2 lw 2
set linetype  3 lc rgb c1 lw 2
set linetype  4 lc rgb c3 lw 2
set linetype  5 lc rgb "black" lw 2
set linetype cycle 4

set style line 1 lt 1 pt 5 ps 0.5 
set style line 2 lt 2 pt 7 ps 0.4 
set style line 3 lt 3 pt 8 ps 0.5 

p="er"
m='plot for [i=1:words(ss)] p.word(ss,i) u "cpr":v w lp ls i t word(ss,i)'
s='set output "a".v.".pdf"' 
set logscale x 2
set xrange [0.5:32]
set xlabel "R / h"

set format y "10^{%L}"
set yrange [1e-30:1]
set logscale y
set ylabel "velocity error L_2"
v="vl2x" ; @s ; @m
set ylabel "We_{max}" offset -1
v="vlmx" ; @s ; @m 

unset format y
unset logscale y
set yrange [0:2]
set ylabel "{/Symbol D}p / p_L"
v="pd" ; @s ; @m

unset yrange
set ylabel "volume"
v="m2" ; @s ; @m

unset format y
unset logscale y
set yrange [-1:1]
set ylabel "center velocity error"
v="v2x" ; @s ; @m

unset format y
unset logscale y
set yrange [-1:1]
set ylabel "center error"
v="c2x" ; @s ; @m
