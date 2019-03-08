#!/usr/bin/env gnuplot

set terminal pdfcairo

f="kerr"


set xrange [3:]
set yrange [0.001:1]
set yrange [0.:0.1]
set xtics 1,2
set output "a.pdf"
#set logscale y
plot \
"kerrcpr4" u 'n':'ekma' w lp pt 7 t "cpr=4",\
"kerrcpr2" u 'n':'ekma' w lp pt 7 t "cpr=2",\
"kerrcpr2a" u 'n':'ekma' w lp pt 7 t "cpr=2a",\
"kerrcpr1" u 'n':'ekma' w lp pt 7 t "cpr=1",\
"kerrcpr1a" u 'n':'ekma' w lp pt 7 t "cpr=1a",\
