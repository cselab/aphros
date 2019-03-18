#!/usr/bin/env gnuplot

set macros
reset

set terminal pdfcairo color font "Helvetica,14" size 3.2,2.2 enhanced

c0='#1f77b4'
c1='#ff7f0e'
c2='#2ca02c'
c3='#d62728'
c4='#9467bd'
c5='#8c564b'

set linetype  1 lc rgb c0 lw 2
set linetype  2 lc rgb c2 lw 2
set linetype  3 lc rgb c1 lw 2
set linetype  4 lc rgb c3 lw 2
set linetype  5 lc rgb c4 lw 2
set linetype  6 lc rgb c5 lw 2
set linetype cycle 6

#set key vertical outside spacing 1.25
unset key

set style line 1 lt 1 pt 7 ps 0.5 pi 5
set style line 2 lt 2 pt 5 ps 0.5 pi 5
set style line 3 lt 3 pt 8 ps 0.5 pi 5
set style line 4 lt 4 pt 2 ps 0.5 pi 5 dt 2
set style line 5 lt 5 pt 2 ps 0.5 pi 5


wall = 0

# frame index to simulation time
sim = 0.01
# bubble radius
r = 0.15

# shift in y
shy = wall ? -r : -0.5
# shift in t
dt = 0.00
dtk = 1/1.06

set key
set key right bottom

# case (directory)
c="center"
ll = system("echo ".c."/*/ch/neck")
tt = system("for f in ".ll." ; do g=${f#".c."/} ; echo ${g%/ch/neck} | tr -d '_' ; done")

set xlabel "t / T"

set ylabel "r_n"
set output "rn.pdf"
krn = 1.
plot \
"ref/rn" w l ls 5  lc "black" t "BI" , \
"ref/rnexp" w p pt 7 ps 0.1 lc "black" t "exp" , \
for [i=1:words(ll)] word(ll,i) u (column('i')*sim*dtk-dt):((column('z1')+shy)/r*krn) w lp ls i t word(tt,i) , \

#for [i=1:words(ll)] word(ll,i) u ($0*sim):(($1+shy)/r) w lp ls i t ""


set ylabel "R_c"
set output "rc.pdf"
set xrange [0:0.3]
set yrange [0:3]

# units from ref/rc*
# capillary time [ms]
tcap = 0.683  # ms
# bubble radius [mm]
rb = 0.3
krc = 1.0

plot \
"ref/rcm" u ($1/tcap):($2/rb) w l lc "black" t "low" , \
"ref/rcp" u ($1/tcap):($2/rb) w l lc "black" t "up" , \
"ref/rcbi" u 1:2 w l lc "red" t "BI" , \
"ref/rcexp" w p pt 7 ps 0.1 lc "black" t "exp" , \
for [i=1:words(ll)] word(ll,i) u (column('i')*sim*dtk-dt):(column("r0")/r*krc) w lp ls i t word(tt,i) , \
for [i=1:words(ll)] word(ll,i) u (column('i')*sim*dtk-dt):(column("r1")/r*krc) w lp ls i t word(tt,i) , \



set xlabel "r_n / R"
set ylabel "R_c / R"
set output "rcrn.pdf"
set xrange [0:1]
set yrange [0:5]
plot \
"ref/rcrn" u 1:2 w l lc "black" t "BI" , \
"ref/rcrnexp" w p pt 7 ps 0.1 lc "black" t "exp" , \
for [i=1:words(ll)] word(ll,i) u ((column("z1")+shy)/r*krn):(column("r0")/r*krc) w lp ls i t word(tt,i) , \

