#!/usr/bin/env gnuplot

set macros
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

if (!exists("yrange")) {
  yrange = "[1e-30:1]"
}

set linetype  1 lc rgb c0 lw 2
set linetype  2 lc rgb c2 lw 2
set linetype  3 lc rgb c1 lw 2
set linetype  4 lc rgb c3 lw 2
set linetype cycle 9

set xlabel "t / T"

p="/traj.dat"
m='plot for [f in ss] f.p u "t":v w l t f'
r=', "ex".p u "t":v w l lc rgb "black" dt 2'
s='set output "a".v.".pdf"'

v="pd"

set ylabel "{/Symbol D}p / p_L"
@s
@m @r

set format y "10^{%L}"
set xrange [0:1]
set yrange @yrange
set logscale y
set ylabel "We_{max}" offset -1

v="vlmx"

@s
wek=system('python -c "\
exec(open('."'".'par.py'."'".').read()) ; \
print(bbr[0][0] * 2 / sig)"')

plot for [f in ss] f.p u "t":(column(v)**2*wek) w l t f
