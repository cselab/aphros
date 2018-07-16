set macros
reset

set terminal pdf color
set output "er.pdf"

set grid

set logscale x 2
set logscale y 2

d="."

set yrange [0.003:1]

plot for [f in "mfer gerris"] \
"<grep ".f." ".d."/err | cut -d' ' -f 2,5" w linespoints pt 5 lw 1 t "".f.""

