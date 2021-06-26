#!/usr/bin/env gnuplot

set terminal qt noenhanced
set macro

V=0.283
f(n) = 2 * pi * (n / 6. - 1)

#p for [f in system('echo sides4_*/')] f.'/traj.dat' u 't':'volume_0' w l lw 2 t f \
#, 'traj.dat' u 't':'volume_0' w l lw 2 t 'current' \
#, for [n in "4"] V+x*f(n) t n ls '--' lw 2 lc 'black' \

sides=5
p for [f in system('echo sides'.sides.'_*/')] '<./calc_rate --sides '.sides.' --traj '.f.'/traj.dat' u 't':'dAdt' w l lw 2 t f \
, '' u 't':'dAdt_exact' w l lw 2 lc 'black'
