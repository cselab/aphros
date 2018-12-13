set terminal pdfcairo

#ss="chnp9 ba ch ge"
ss="chnp9 chnp3 ch"
reset
set output "aem.pdf"
set logscale y
plot for [f in ss] f."/kerr" u 'cpr':'em' w p t f

reset
set output "ae2.pdf"
set logscale y
plot for [f in ss] f."/kerr" u 'cpr':'e2' w p t f
