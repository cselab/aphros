#!/usr/bin/env gnuplot

set macro

set terminal pdf noenhanced

ll1="\
pasc_nn096 \
pasc_nn192 \
pasc_nn384 \
"

ll2="\
pasc_nn096_vof \
pasc_nn192_vof \
pasc_nn384_vof \
"

ll=ll2

nn="200 200 200"

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

hist=sprintf("ch.hist --range %g %g --bins 100", hx2*0.25, hx2*128)
p for [i=1:words(ll)] \
  "<o=hist_".word(ll,i)." ; ( test -f $o && cat $o ) || (./trajcol ".word(ll,i)."/traj_0{170..199}.csv | ".hist." ) | tee $o"  \
  u ($1/hx2):($2/30) \
  w l lw 2 t word(ll,i) ,\
  (x/32)**(-10/3) lc "black"

