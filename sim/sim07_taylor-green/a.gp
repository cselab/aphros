set macros
reset

#set terminal pdf dashed color
set terminal pdf color
set output "kedr.pdf"

set grid

v='dk'

l="nx64 nx128p64 nx128p512t20 nx128p64t20i2"

plot for [f in l] "<paste ".f."/sc/{t,".v."}" w l lw 3 t "".f.""

set xrange [0:2]
set yrange [-0.002:0.014]

