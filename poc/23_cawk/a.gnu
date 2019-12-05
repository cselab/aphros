#!/usr/bin/env gnuplot

set macros
set terminal pdfcairo color enhanced

set output "a.pdf"

#p  \
#"stat" u "r0":"r1" w p ps 0.6 pt 7 , x

#"<awk 'NR==1 || $8 > 0.0000' stat" u "t":(column("v0")/column("v")) w p ps 0.6 pt 7 ,\
#"<awk 'NR==1 || $8 > 0.0000' stat" u "t":(column("v1")/column("v")) w p ps 0.6 pt 7 ,\

#exit

#set size 1
#p [:2][:1] "<awk 'NR==1 || $5 > 0.0001' stat" u "x":"y":((column("v1"))) w circle palette

#set style fill  transparent solid 0.35 noborder
#set style circle radius 0.01

#set pal neg gray

#set cbrange [-14:0]
#p [:2][:1] "<awk 'NR==1 || $5 > 0.0001' stat" u "x":"y":((column("v1"))) w circle palette
#exit

p "<awk 'NR>1 {print $8}' stat | ch.hist" w l

#p  \
#"<./a.awk -v a=8 stat " smooth cspline , \
