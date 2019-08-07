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
  ss = "ch ba"
}

set linetype  1 lc rgb c0 lw 2
set linetype  2 lc rgb c2 lw 2
set linetype  3 lc rgb c1 lw 2
set linetype  4 lc rgb c3 lw 2
set linetype cycle 9

set style line 1 lt 1 pt 7 ps 0.5
set style line 2 lt 2 pt 5 ps 0.5
set style line 3 lt 3 pt 8 ps 0.5 

set xlabel "t / T"

p="/traj.dat"
m='plot for [f in ss] f.p u "t":v w l t f'
r=', "ex".p u "t":v w l lc rgb "black" dt 2 t "flow"'
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

set format y "10^{%L}"
set xrange [0:1]
set yrange [1e-7:1]
set logscale y
set ylabel "We_{max}" offset -1

v="vlmx"

@s
wek=system('python -c "\
exec(open('."'".'par.py'."'".').read()) ; \
print(bbr[0][0] * 2 / sig)"')

plot for [i=1:words(ss)] word(ss,i).p u "t":(column(v)**2*wek) w l ls i t word(ss,i)

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

