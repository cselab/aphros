#!/usr/bin/env gnuplot

set terminal pdfcairo
set output "rate.pdf"
set macro

V=0.283
f(n) = 2 * pi * (n / 6. - 1)

sides=4
p for [f in system('echo sides'.sides.'_*/')] '<./calc_rate --sides '.sides.' --traj '.f.'/traj.dat' u 't':'dAdt' w l lw 2 t f \
, '' u 't':'dAdt_exact' w l lw 2 lc 'black'
