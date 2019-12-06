#!/usr/bin/env gnuplot

set macro

set terminal pdf noenhanced

ll="\
pasc_nn096 \
pasc_nn192 \
pasc_nn384 \
"

nn="160 160 160 200"

hx0=1./96
hx1=1./192
hx2=1./384
set xrange [0.5:32]
set yrange [0.1:10000]

d="/mnt/daint/scratch/snx3000/karnakov/fountain/"

set xlabel "radius / finest mesh cell

v="hist"
set log x 2
set log y 10
set output v.".pdf"
set arrow from hx0,graph 0 to hx0,graph 1 nohead
set arrow from hx1,graph 0 to hx1,graph 1 nohead
set arrow from hx2,graph 0 to hx2,graph 1 nohead

hist=sprintf("ch.hist --range %g %g --bins 100", hx2*0.25, hx2*32)
p for [i=1:words(ll)] \
  "<o=hist_".word(ll,i)." ; ( test -f $o && cat $o ) || (cd ".d." && ch.gettrajcol r ".word(ll,i)."/traj_0".word(nn,i).".csv | ".hist." ) | tee $o"  \
  u ($1/hx2):2 \
  w l lw 2 t word(ll,i) ,\
  (x/30)**(-10/3) lc "black"


exit

#!  ./radius daintthickvof/traj_01??.csv  > qvof

#! ./radius @d/waterfall_vofm_Nov24_wip_thick/traj_020?.csv  > q
! ./radius @d/waterfall_vofm_Nov25_wip_u1/traj_011?.csv  > q


reset
set log
set output "h2.pdf"
set xrange [0.01:0.1]
set yrange [0.1:1000]
p \
  "< <q ch.hist"  w l ,\
  "< <qvof ch.hist"  w l ,\
  (x/0.045)**(-10./3) lc "black" 
exit



