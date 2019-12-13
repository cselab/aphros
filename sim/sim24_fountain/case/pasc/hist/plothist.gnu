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

ll=ll1.ll2

hx2=1./384
set yrange [0.01:1000]

d="/mnt/daint/scratch/snx3000/karnakov/fountain/"

v="hist"
set log x 2
set log y 10
set output v.".pdf"

physx = 100  # factor: mm / simulation length


binL = 0
binR = 0.1
tt = "traj_0{201..240}.csv"
bins = 50
files = int(system("echo ".tt." | wc -w"))
hist=sprintf("ch.hist --range %g %g --bins %d | awk '{print $1 , $2/%d}'", binL, binR, bins, files)
info = sprintf('bin width(mm)=%g, bins=%d, files=%d', (binR - binL)/bins*physx, bins, files)
system("echo '".info."' > bins")

set xrange [binL*100:binR*100]
p for [i=1:words(ll)] \
  "<o=hist_".word(ll,i)." ; ( test -f $o ) || (./trajcol ".word(ll,i)."/".tt." | ".hist." > tmp ; mv -f tmp $o ) && cat $o"  \
  u ($1*100):($2/30) \
  w lp pt 7 ps 0.5 lw 2 t word(ll,i) ,\
  (x/1.6)**(-10/3) lc "black"

