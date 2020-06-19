set terminal pdf

set output "a.pdf"
p "<ap.gettraj 0" u "n":"x" w lp
