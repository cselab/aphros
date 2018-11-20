set macros
reset

set terminal pdfcairo color font ",18" size 5,2.5 enhanced

set nokey
c0='#1f77b4'
c1='#ff7f0e'
c2='#2ca02c'
c3='#d62728'
c4='#9467bd'

set linetype  1 lc rgb c0 lw 3 
set linetype  2 lc rgb c2 lw 3
set linetype  3 lc rgb c1 lw 3
set linetype  4 lc rgb c3 lw 3
set linetype cycle 9

p="/traj.dat"
ll="ch ba"
m='plot for [f in "ch ba"] f.p u'

set output "ax.pdf"
@m "t":"c2x" w l

set output "ay.pdf"
@m "t":"c2y" w l

set output "az.pdf"
@m "t":"c2z" w l

set output "ap.pdf"
@m "t":"pd" w l

set output "avlmx.pdf"
@m "t":"vlmx" w l

set output "avlmy.pdf"
@m "t":"vlmy" w l

set output "avlmz.pdf"
@m "t":"vlmz" w l

set output "avl2x.pdf"
@m "t":"vl2x" w l

set output "avl2y.pdf"
@m "t":"vl2y" w l

set output "avl2z.pdf"
@m "t":"vl2z" w l

