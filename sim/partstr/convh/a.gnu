#!/usr/bin/env gnuplot

set terminal pdfcairo
set macro

f="kerr"


set xrange [2:]
#set yrange [0.001:1]
#set yrange [0.:0.3]
set output "a.pdf"
set logscale y
t='(1. + column("n") * 0.25)'
plot \
f u @t:'ek2' w lp pt 7 t 'exact,L2' ,\
f u @t:'ekm' w lp pt 7 t 'exact,max' ,\
#f u @t:'ek2a' w lp pt 7 t 'exactavg,L2' ,\
#f u @t:'ekma' w lp pt 7 t 'exactavg,max' ,\
