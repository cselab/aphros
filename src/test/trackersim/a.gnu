#!/usr/bin/env gnuplot

set terminal pdf

set output "a.pdf"
p \
"<ap.gettraj 0" u "n":"x" w lp , \
"<(cd '".ARG1."' && ap.gettraj 0)" u "n":"x" w lp t "ref"
