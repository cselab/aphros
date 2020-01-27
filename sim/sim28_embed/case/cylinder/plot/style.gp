fnt="Times,14"
if (!exists("termsize")) termsize="2.15,2"
set terminal pdfcairo color font fnt size @termsize enhanced
lbx = "{/:Italic x}"
lbvx = "{/:Italic u_{/:Italic x}}"
lbt = "{/:Italic t}"

set tics scale 0.5
set xtics offset 0,0

set xlabel offset 0,0.25

set lmargin 8
set rmargin 2
set tmargin 0.5
set bmargin 3.5

set macro
set nokey

c0='#1f77b4'
c1='#ff7f0e'
c2='#2ca02c'
c3='#d62728'
c4='#9467bd'
c9='#000000'

set linetype  1 lc rgb c0 lw 2
set linetype  2 lc rgb c1 lw 2
set linetype  3 lc rgb c2 lw 2
set linetype  4 lc rgb c9 lw 2
set linetype  5 lc rgb c4 lw 2
set linetype cycle 5

set dashtype 2 (11,6)
set dashtype 3 (12,6,4,6)

set style line 1 lt 1 pt 7 ps 0.5
set style line 2 lt 2 pt 5 ps 0.5 dt 2
set style line 3 lt 3 pt 9 ps 0.5 dt 3
set style line 4 lt 4 pt 2 ps 0.5

