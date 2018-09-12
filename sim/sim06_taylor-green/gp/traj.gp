set macros
reset

#set terminal pdf dashed color
set terminal pdf color
set output "traj.pdf"

set grid

v='x'

l=system("ls nx128*.dat")

set xlabel "n"
set ylabel v
plot for [f in l] "<awk 'NR!=2' ".f using "n":"".v."" w l lw 3 t "".f.""

set xrange [0:2]
set yrange [-0.002:0.014]

