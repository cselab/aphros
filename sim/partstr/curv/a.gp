set macros
reset

#set terminal pdf dashed color
set terminal pdf color
set output "er.pdf"

set grid

set logscale x 2
set logscale y 10

plot for [f in "mfer gerris"] \
"<grep ".f." err | cut -d' ' -f 2,5" w linespoints pt 5 lw 1 t "".f.""
