#!/usr/bin/env gnuplot

set terminal pdf

set output "a.pdf"
p \
"c2x" w lp t "current" , \
ARG1."/c2x" w lp t "ref"
