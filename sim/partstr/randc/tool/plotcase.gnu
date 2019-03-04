#!/usr/bin/env gnuplot

set terminal pdfcairo color font "Helvetica,14" size 3.2,2.2 enhanced
set macro
reset

set key font ",10"

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
set linetype cycle 5

set style line 1 lt 1 pt 7 ps 0.5
set style line 2 lt 2 pt 5 ps 0.5
set style line 3 lt 3 pt 8 ps 0.5
set style line 4 lt 4 pt 2 ps 0.5

c='ls -d dim3_overlap1_segcirc1_np9_ns?_hp4'
c='ls -d dim3_overlap1_segcirc1_np9_ns?_hp4'
c='ls -d dim3_overlap1_segcirc1_np9_ns2_hp?'
c='ls -d dim3_overlap1_segcirc?_np9_ns2_hp4'
c='ls -d dim2_overlap0_segcirc?_np9_ns2_hp4'

if (!exists("ss")) {
  ss = system("ls -d dim*")
}

if (!exists("pre")) {
  pre = "a"
}

tl='system("tr -d _ <<< ".word(ss,i))'
s='pt 7 ps 0.5'
m='for [i=1:words(ss)] word(ss,i)."/ch/kerravg" \
    u "cpr":@v w lp t '.tl.' ls i'
mfs='for [i=1:words(ss)] word(ss,i)."/ch/kerravg" \
    u "cpr":@v."sl":@v."sh" w filledcurve t "" ls i fs transparent solid 0.5'
mf='for [i=1:words(ss)] word(ss,i)."/ch/kerravg" \
    u "cpr":@v."l":@v."h" w filledcurve t "" ls i fs transparent solid 0.25'

o='set output pre.@v.".pdf"'

set logscale x 2

set xrange [0.5:32]
set xlabel "R / h"

set logscale y
set ylabel "L_{/Symbol \245}"
set yrange [0.001:10]
v='"em"'
@o ; plot @m , @mfs

set logscale y
set ylabel "L_2"
set yrange [0.001:10]
v='"e2"'
@o ; plot @m  , @mfs

set logscale y
set ylabel "L_s"
set yrange [0.001:10]
v='"es"'
@o ; plot @m  , @mfs


unset logscale y
set yrange [0:3]
set ylabel "V / V_{exact}"
v='"vf"'
@o ; plot @m  , @mfs
