#!/usr/bin/env gnuplot
  #

load "style.gnu"

d="stat.dat"

set output "a.pdf"

set xlabel "t / T"
set ylabel "perimeter fluctiation"

nx = 128
r = 0.1

set xtics 100
unset key

# initial perimeter
per0 = 2 * pi * r

# equivalent radius
Freq = '((column("m2")*nx/pi)**0.5)'
# perimeter
Fper = '(column("area")*nx)'

p d u "t":((@Fper - @Freq * 2*pi)/per0) every ::1 w l lw 0.2

sig = 0.000164493406685
rho1 = 1
freq = 1

# bubble radius for which mode n oscillates with frequency freq
R(n) = (n * (n - 1) * (n + 1) * sig / (rho1 * (2 * pi * freq * 0.5) ** 2)) ** (1./3)

set output "b.pdf"
reset
set xlabel "R_{eq} / R_0"
set ylabel "perimeter fluctiation"
unset key

vl(x) = sprintf('set arrow from %g, graph 0 to %g, graph 1 nohead front', x, x)

do for [n=4:8] {
  #x=R(n)*0.96
  x=R(n) / r
  eval(vl(x))
  set label sprintf("%d", n) at x,0.2 offset 0.4
}

p d u (@Freq / r):((@Fper - @Freq * 2*pi)/per0) every ::1 w l lw 0.2


set output "c.pdf"
reset
set xlabel "t / T"
set ylabel "(V - V_0) / V_0"
unset key

p d u "t":"m2d" every ::1 w l lw 0.2

