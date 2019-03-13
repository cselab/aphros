#!/usr/bin/env gnuplot

set macros
reset

set terminal pdfcairo color font "Helvetica,14" size 3.2,2.2 enhanced

c0='#1f77b4'
c1='#ff7f0e'
c2='#2ca02c'
c3='#d62728'
c4='#9467bd'
c5='#8c564b'

set linetype  1 lc rgb c0 lw 2
set linetype  2 lc rgb c2 lw 2
set linetype  3 lc rgb c1 lw 2
set linetype  4 lc rgb c3 lw 2
set linetype  5 lc rgb c4 lw 2
set linetype  6 lc rgb c5 lw 2
set linetype cycle 6

#set key vertical outside spacing 1.25
unset key

set style line 1 lt 1 pt 7 ps 0.5 pi 10
set style line 2 lt 2 pt 5 ps 0.5 pi 10
set style line 3 lt 3 pt 8 ps 0.5 pi 10
set style line 4 lt 4 pt 2 ps 0.5 pi 10
set style line 5 lt 5 pt 2 ps 0.5 pi 10

set style line 11 lt 1 pt 6 ps 0.5 pi 10 dt 2
set style line 12 lt 2 pt 4 ps 0.5 pi 10 dt 2
set style line 13 lt 3 pt 8 ps 0.5 pi 10 dt 2


# frame index to simulation time
sim = 0.01
# bubble radius
r = 0.15

unset key

set xrange [0:0.4]

set xtics 0.1
set ytics 0.2

#ll = system("ls -d *symm1* | grep -v rho3k")
#ll = system("echo nx512*symm1*")
llc = system("echo nx{064,128,256}*symm1/ch/neck")
llb = system("echo nx{064,128,256}*symm1/ba/neck")
tt = "64 128 256"
# shift in y
shy = -0.5


m='plot \
"ref/rnexp" w p pt 7 ps 0.1 lc "black" t "exp" , \
"ref/rn" w l lc "black" t "BI" , \
  for [i=1:words(ll)] word(ll,i) u ($0*sim):(($2+shy)/r) w lp ls i t word(tt,i)'

set xlabel "t / T"
set ylabel "r_n"

if (exists("wall")) {
  llc = system("echo wall/*/ch/neck")
  tt = system("for f in ".llc." ; do g=${f%ns*ch/*} ; echo ${g#wall/nx} ; done")
  set key
  set xrange [0:5]
  set xtics auto
  set ytics auto
  set style line 1 lt 1 pt 7 ps 0.5 pi 50
  set style line 2 lt 2 pt 5 ps 0.5 pi 50
  set style line 3 lt 3 pt 8 ps 0.5 pi 50
  shy = -r
  m='plot \
  for [i=1:words(ll)] word(ll,i) u ($0*sim):(($2+shy)/r) w lp ls i t word(tt,i) , \
  for [i=1:words(ll)] word(ll,i) u ($0*sim):(($1+shy)/r) w lp ls i t "" , \
    '
}

set output "rnch.pdf"
ll = llc
@m ;

set output "rnba.pdf"
ll = llb
@m ;

exit

# divided by t**0.5
set output "b.pdf"
set yrange [1:2]
plot \
  "ref/rnexp" u 1:($2/$1**0.5) w p pt 7 ps 0.1 lc "black" t 'exp' , \
  "ref/rn" u 1:($2/$1**0.5) w l lc "black" t 'BI' , \
  for [i=1:words(ll)] word(ll,i).sf u ($0*sim):(($2-0.5)/r/($0*sim)**0.5) w lp ls i t word(tt,i)

