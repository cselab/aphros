#!/usr/bin/env gnuplot

#!rm -f {vx,vy}_*.{xmf,h5}
#!./rs

termsize="3,2"
load "style.gp"

#St = 0.177
St = 0.178

print(sprintf("St=0.164   CD=1.34 ±0.011   CL=±0.315"))

set key font ",8"

flt(x) = system("echo '".x."' | tr '_' '-'")

set key bottom

set xlabel lbt
set xrange [183:200]

nn="-E '26'"

set fit quiet
set fit logfile '/dev/null'

dat="cylinder_re100_ny512/stat.dat"
dat="cylinder_re100_ext40_ny512/stat.dat"

f(x) = sin(4*pi*(x-t0)*St) * Ax + CD
fit [180:200] f(x) dat u "t":(column("dragx")*2) via CD,St,t0,Ax

set output "drag.pdf"
set ylabel "{/:Italic C}_D" offset 0.8
set ytics 0.01
plot for [f in dat]\
  f u "t":(column("dragx")*2) w l t flt(f) ,\
  f(x)

g(x) = sin(4*pi*(x-t0-t0y)*St*0.5)*Ay
fit [180:200] g(x) dat u "t":(column("dragy")*2) via t0y,Ay

print(dat)
print(sprintf("St=%.3f   CD=%.2f ±%.3f   CL=±%.3f", St, CD, Ax, Ay))

set output "lift.pdf"
set ylabel "{/:Italic C}_L"
set ytics 0.1
plot for [f in dat]\
  f u "t":(column("dragy")*2) w l t flt(f) ,\
  g(x)

  #sin(2*pi*(x-181.49)*St)*0.43
