fnt="Helvetica,14"
if (!exists("termsize")) termsize="3.2,2.2"
set terminal pdfcairo color font fnt size @termsize enhanced
lbx = "{/:Italic x}"
lbC = "{/:Italic C}"
lby = "{/:Italic y}"
lbt = "{/:Italic t}"
lw = 1.5

set lmargin 10
set rmargin 12

set macro
set nokey

c0='#1f77b4'
c1='#ff7f0e'
c2='#2ca02c'
c3='#d62728'
c4='#9467bd'
c5='#8c564b'
c9='#000000'

set linetype  1 lc rgb c0 lw lw
set linetype  2 lc rgb c2 lw lw
set linetype  3 lc rgb c1 lw lw
set linetype  4 lc rgb c3 lw lw
set linetype  5 lc rgb c4 lw lw
set linetype  6 lc rgb c5 lw lw
set linetype  7 lc rgb c9 lw lw
set linetype cycle 6

set style line 1 lt 1 pt 7 ps 0.5 pi 100
set style line 2 lt 2 pt 5 ps 0.5 pi 100
set style line 3 lt 3 pt 8 ps 0.5 pi 100
set style line 4 lt 4 pt 2 ps 0.5 pi 100
set style line 5 lt 5 pt 2 ps 0.5 pi 100
set style line 6 lt 6 pt 2 ps 0.5 pi 100
set style line 7 lt 7 pt 2 ps 0.5 pi 100
