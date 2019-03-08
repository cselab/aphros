#!/usr/bin/env gnuplot

set terminal pdfcairo

f="kerr"


set xrange [1:]
set yrange [0.001:1]
set xtics 1
set output "a.pdf"
set logscale y
plot \
f u 'n':'er2' w lp pt 7 t 'finest,L2' ,\
f u 'n':'ek2' w lp pt 7 t 'exact,L2' ,\
f u 'n':'erm' w lp pt 7 t 'finest,max' ,\
f u 'n':'ekm' w lp pt 7 t 'exact,max' ,\
f u 'n':'std' w lp pt 7 t 'std' ,\
