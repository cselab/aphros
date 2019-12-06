#!/usr/bin/env gnuplot

set macros
set terminal pdfcairo color enhanced

set output "a.pdf"


set datafile separator ","
a='a_0158.csv'
b='a_0159.csv'

#h=1./192
#p a u "cl":"r" w p pt 7 ps 0.1 , b u "cl":"r" w p pt 7 ps 0.1
p \
"<./match ".a." ".b u "x":"r" w p pt 7 ps 0.1 ,\
"<./match ".b." ".a u "x":"r" w p pt 7 ps 0.1
