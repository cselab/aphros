set macros
reset

#set terminal pdf dashed color
set terminal pdf color
set output "kedr.pdf"

set grid

v='dk'

l=system("ls -d nx128*/sc/")

plot for [f in l] "<paste ".f."/{t,".v."}" w l lw 3 t "".f.""

set xrange [0:2]
set yrange [-0.002:0.014]

