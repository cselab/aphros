#!/usr/bin/env gnuplot

set terminal pdfcairo

f="kerr"

unset key

set output "aer.pdf"
plot \
f u 'eta':'ek2' w lp pt 7 t 'error' ,\


set output "ait.pdf"
set yrange [0:]
plot \
f u 'eta':'it' w lp pt 7 t 'iter' ,\
"" u 'eta':'it0':'it1' with filledcurve fs transparent solid 0.5 t '' lt 1 ,  \
