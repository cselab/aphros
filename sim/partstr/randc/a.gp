set terminal pdfcairo color font "Helvetica,14" size 3.2,2.2 enhanced
set macro
reset

c0='#1f77b4'
c1='#ff7f0e'
c2='#2ca02c'
c3='#d62728'
c4='#9467bd'

set linetype  1 lc rgb c0 lw 2
set linetype  2 lc rgb c2 lw 2
set linetype  3 lc rgb c1 lw 2
set linetype  4 lc rgb c3 lw 2
set linetype  5 lc rgb c4 lw 2
set linetype cycle 5

set style line 1 lt 1 pt 7 ps 0.5
set style line 2 lt 2 pt 5 ps 0.5
set style line 3 lt 3 pt 8 ps 0.5

ss="ch ba ge"
#ss="chnp9 chnp3 chol"

s='pt 7 ps 0.5'
m='for [i=1:words(ss)] word(ss,i)."/kerravg" \
    u "cpr":@v w lp t word(ss,i) ls i'
mfs='for [i=1:words(ss)] word(ss,i)."/kerravg" \
    u "cpr":@v."sl":@v."sh" w filledcurve t "" ls i fs transparent solid 0.5'
mf='for [i=1:words(ss)] word(ss,i)."/kerravg" \
    u "cpr":@v."l":@v."h" w filledcurve t "" ls i fs transparent solid 0.25'

set logscale x 2

set xrange [0.5:16]
set xlabel "cpr"

set output "aem.pdf"
set logscale y
set ylabel "L_{/Symbol \245}"
set yrange [0.01:10]
v='"em"'
plot @m , @mfs , @mf

set output "ae2.pdf"
set logscale y
set ylabel "L_2"
set yrange [0.01:10]
v='"em"'
v='"e2"'
plot @m , @mfs , @mf

set output "avf.pdf"
unset logscale y
set yrange [0:3]
set ylabel "V / V_{exact}"
v='"vf"'
plot @m , @mfs , @mf
