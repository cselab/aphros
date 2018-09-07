set macros
reset

set terminal pdf color
set output "a.pdf"

set grid

v='x'

l="r1 r10 r100 ''"

set xlabel "n"
set ylabel v
plot for [f in l] "traj".f."/dat/trajvx2_mfer.dat" w l lw 3 t f

set xrange [0:2]
set yrange [-0.002:0.014]

