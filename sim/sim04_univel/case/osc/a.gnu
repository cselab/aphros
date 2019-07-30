#!/usr/bin/env gnuplot

set macros
reset

set terminal pdfcairo color font "Helvetica,14" size 3.2,2.2 enhanced
set terminal pdfcairo color font "Helvetica,14" size 4.2,2.2 enhanced

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

set style line 1 lt 1 pt 7 ps 0.5 pi 5
set style line 2 lt 2 pt 5 ps 0.5 pi 5
set style line 3 lt 3 pt 8 ps 0.5 pi 5
set style line 4 lt 4 pt 2 ps 0.5 pi 5 dt 2
set style line 5 lt 5 pt 2 ps 0.5 pi 5

# frame index to simulation time
sim = 10./100
# bubble radius
r = 0.3

set key outside

ll = system("echo */{ch,ba}")
tt = system("for f in ".ll." ; do echo $f | tr -d '_' ; done")

# units from simulation
rb = 0.3
tb = 1
# mode n=2 period
n = 2
R = 1.1 ** (1. / 3.)
tn = 2 * pi / ((n-1)*(n+1)*(n+2)) ** 0.5 * (R**3) ** 0.5
s2 = n * (n + 1) * (n - 1) * (n + 2) / ((n+1)*0.010 + n * 1)

tn = 2 * pi / s2 ** 0.5
#tn = 2 * pi / ((n-1)*(n+1)*(n+2)) ** 0.5 * (R**3) ** 0.5


la = "nx064_ns2_symm0_wall0_tmax10/"
lb = "nx128x_ns2_symm0_wall0_tmax10/"
l3 = "nx256_ns2_symm0_wall0_tmax10/"
l4 = "nx256_ns2_symm0_wall0_tmax10_r075/"

set output "rz.pdf"
set xlabel "t / t_b"
set ylabel "{/Symbol D} r_z"
unset xrange
unset yrange

tnlk = 1.07
#tnl = tn * 1.12

set samples 1000

plot \
la."ch/neck" u ((column('i')*sim)):((0.5-column("z0"))/rb-1) w l t "r0.30-ch"  \
, la."ba/neck" u ((column('i')*sim)):((0.5-column("z0"))/rb-1) w l t "r0.30-ba"  \
, lb."ch/neck" u ((column('i')*sim)):((0.5-column("z0"))/rb*2.-1) w l t "r0.15-ch"  \
, lb."ba/neck" u ((column('i')*sim)):((0.5-column("z0"))/rb*2.-1) w l t "r0.15-ba"  \
, exp(-x*0.1)*cos(pi+2*pi*((x-0.4))/(tn*tnlk))*0.03+0.032 w l lc rgb "black" lw 1 t sprintf("lin,T*%g", tnlk)  \

set output "ekin.pdf"
set ylabel "E_{kin}"

plot \
la."ch/stat.dat" u (column('t')):(column("ekin")*1.57e3) w l t "r0.30-ch" \
, lb."ch/stat.dat" u (column('t')):(column("ekin")*5e4) w l t "r0.15-ch" \
, l3."ch/stat.dat" u (column('t')):(column("ekin")*2e5) w l t "r0.15x-ch" \
, l4."ch/stat.dat" u (column('t')):(column("ekin")*6.25e6) w l t "r0.075-ch" \
, exp(-x*0.36)*(sin(pi+2*pi*((x-0.0))/(tn*tnlk))*0.287)**2 w l lc rgb "black" lw 1 t  sprintf("lin,T*%g", tnlk)  \


