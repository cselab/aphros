#!/usr/bin/env gnuplot

V=0.283
f(n) = 2 * pi * (n / 6. - 1) * 2

p 'traj.dat' u 't':'volume_0' w l lw 2 \
, for [n in "2 3 4 6 8"] V + x * f(n) t n \
