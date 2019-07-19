#!/usr/bin/env gnuplot

  #

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

set style line 1 lt 1 pt 7 ps 0.5 pi 5
set style line 2 lt 2 pt 5 ps 0.5 pi 5
set style line 3 lt 3 pt 8 ps 0.5 pi 5
set style line 4 lt 4 pt 2 ps 0.5 pi 5 dt 2
set style line 5 lt 5 pt 2 ps 0.5 pi 5


# case (directory)
c="wall" ; wall = 1
#c="center" ; wall = 0

# frame index to simulation time
sim = 0.01
# bubble radius
r = 0.15
r2 = 0.075


# shift in y
shy = wall ? -r : -0.5
shy2 = wall ? -r2 : -0.5
# shift in t
dt = 0
dtk = 1.2
#dtk = 1.14
#dtk = 1

set key
set key right bottom

ll = system("echo ".c."/nx*256*/ch/neck")
tt = system("for f in ".ll." ; do g=${f#".c."/} ; echo $g | sed -r 's/(nx...).*_(:?r075)?.*\\/ch.*/\\1\\2/' | tr -d '_' ; done")

l1 = word(ll, 1)
l2 = word(ll, 2)
l3 = word(ll, 3)

ii1 = '1:2'
ii2 = '3:3'

set xlabel "t / T"
set xrange [0.:0.4]
#unset yrange
#set logscale x 2
#set logscale y 2

set output "rn.pdf"
set xtics 0.1
set ylabel "r_n / R"
set yrange [0.:1.]
krn = 1.
plot \
"ref/rnexp" w p pt 7 ps 0.1 lc "black" t "exp" , \
for [i=@ii1] word(ll,i) u (column('i')*sim*dtk-dt):((column('z1')+shy)/r*krn) w l ls i t word(tt,i) , \
for [i=@ii2] word(ll,i) u (column('i')*sim*dtk-dt):((column('z1')+shy2)/r2*krn) w l ls i t word(tt,i) , \


#"ref/rn" w l ls 5  lc "black" t "BI" , \

set ylabel "R_c / R"
set output "rc.pdf"
set xrange [0:0.3]
set yrange [0:3]

# units from experiment, ref/rc*
# capillary time [ms]
etcap = 0.683  # ms
# bubble radius [mm]
er = 0.3
krc = 1.2

plot \
"ref/rcm" u ($1/etcap):($2/er) w l lc "black" t "low" , \
"ref/rcp" u ($1/etcap):($2/er) w l lc "black" t "up" , \
"ref/rcexp" w p pt 7 ps 0.1 lc "black" t "exp" , \
for [i=@ii1] word(ll,i) u (column('i')*sim*dtk-dt):(column("r0")/r*krc) w l ls i t word(tt,i) , \
for [i=@ii1] word(ll,i) u (column('i')*sim*dtk-dt):(column("r1")/r*krc) w l ls i t word(tt,i) , \
for [i=@ii2] word(ll,i) u (column('i')*sim*dtk-dt):(column("r0")/r2*krc) w l ls i t word(tt,i) , \

#"ref/rcbi" u 1:2 w l lc "red" t "BI" , \

set xlabel "r_n / R"
set ylabel "R_c / R"
set output "rcrn.pdf"
set xrange [0:1]
set yrange [0:5]
plot \
"ref/rcrn" u 1:2 w l lc "black" t "BI" , \
"ref/rcrnexp" w p pt 7 ps 0.1 lc "black" t "exp" , \
for [i=@ii1] word(ll,i) u ((column("z1")+shy)/r*krn):(column("r0")/r*krc) w l ls i t word(tt,i) , \
for [i=@ii2] word(ll,i) u ((column("z1")+shy2)/r2*krn):(column("r0")/r2*krc) w l ls i t word(tt,i) , \


# units from simulation
rb = r * 2 ** (1. / 3)
rb2 = r2 * 2 ** (1. / 3)
tb = 1 * 2 ** 0.5
# mode n=2 period
n = 2
tn = 2 * pi / ((n-1)*(n+1)*(n+2)) ** 0.5

set output "rth0.pdf"
set xlabel "t / t_b"
set ylabel "r_{{/Symbol q}=0} / R_b"
unset xrange
unset yrange
set xtics auto
set ytics 0.2
set key top
set xrange [0:6]
unset yrange
plot \
for [i=@ii1] word(ll,i) u ((column('i')*sim-dt)/tb):((0.5-column("x0"))/rb) w l ls i t word(tt,i) \
, for [i=@ii2] word(ll,i) u ((column('i')*sim-dt)/tb):((0.5-column("x0"))/rb2) w l ls i t word(tt,i) \
, "ref/rth0" u ($1*dtk):2 w l lc "black" t "exp,k=1.2" \
, "ref/rth0" u ($1+0.5):2 w l lc "black" dt '-' t "exp,k=1,shift" \

#cos(pi*0.5+2*pi*((x-2.75)/1.13)/tn)*0.2+1 w l dt 2 lc rgb c4 lw 1 t "lin,T*1.13" \
#, cos(pi*0.5+2*pi*((x-2.25)/0.93)/tn)*0.2+1 w l lc rgb c2 lw 1 t "lin,T*0.93" \
#, cos(pi*0.5+2*pi*((x-2.25))/tn)*0.2+1 w l lc rgb c4 lw 1 t "lin" \
