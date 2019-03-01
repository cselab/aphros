#!/usr/bin/env gnuplot

set macros
reset

set terminal pdfcairo color font "Helvetica,14" size 3.2,2.2 enhanced

set nokey
c0='#1f77b4'
c1='#ff7f0e'
c2='#2ca02c'
c3='#d62728'
c4='#9467bd'

set linetype  1 lc rgb c0 lw 3 
set linetype  2 lc rgb c2 lw 3
set linetype  3 lc rgb c1 lw 3
set linetype  4 lc rgb c3 lw 3
set linetype cycle 9


if (!exists("ll")) {
  ll = "ch ba"
}

p="/traj.dat"
m='plot for [f in ll] f.p u "t":v w l'
r=', "ex".p u "t":v w l lc rgb "black" dt 2'
s='set output "a".v.".pdf"'

v="c2x"
@s
@m @r

v="c2y"
@s
@m @r

v="c2z"
@s
@m @r

v="pd"
@s
@m @r

set logscale y
set yrange [1e-4:10]

v="vlmx"
@s
@m

v="vlmy"
@s
@m

v="vlmz"
@s
@m

v="vl2x"
@s
@m

v="vl2y"
@s
@m

v="vl2z"
@s
@m

