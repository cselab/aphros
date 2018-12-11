set terminal pdfcairo

f="kerr"

reset
set output "a.pdf"
set logscale y
plot f u 1:2 w lp pt 7 t 'ref' , f u 1:3 w lp pt 7 t 'exact'
