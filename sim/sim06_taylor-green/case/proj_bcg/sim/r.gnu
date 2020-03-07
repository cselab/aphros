#!/usr/bin/env gnuplot

termsize = "3,2"
load "style.gp"

set xr[0:20]
set ytics 0.005
set xlabel lbt
set ylabel "energy dissipation rate"

set key outside

ll = system("echo nx*")
tt = "128 256 512"

set output "kedr.pdf"
plot \
   "<paste ../ref/{t,dk}" w l t "ref  " ls 7 , \
   for [i=1:words(ll)] "<paste ".word(ll, i)."/sc/{t,dk}" w l ls i dt 1 t word(tt, i) , \
