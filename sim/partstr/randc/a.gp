set terminal pdfcairo

f="kerr"

reset
set output "a.pdf"
set logscale y
plot f u 3:4 w p 
