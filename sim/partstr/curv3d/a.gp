set macros
reset

set terminal pdf color

set grid

set logscale x 2
set logscale y 2

d="."

set yrange [0.003:1]

set output "ermax.pdf"
plot for [f in "mfer gerris"] \
"<grep ".f." ".d."/err | cut -d' ' -f 2,4" w linespoints pt 5 lw 1 t "".f.""

set output "erl1.pdf"
plot for [f in "mfer gerris"] \
"<grep ".f." ".d."/err | cut -d' ' -f 2,5" w linespoints pt 5 lw 1 t "".f.""

set output "erl2.pdf"
plot for [f in "mfer gerris"] \
"<grep ".f." ".d."/err | cut -d' ' -f 2,6" w linespoints pt 5 lw 1 t "".f.""

