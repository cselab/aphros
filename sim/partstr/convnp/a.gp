set terminal pdfcairo
set output "a.pdf"
set logscale y
plot "kerr" u 1:4 w lp pt 7
