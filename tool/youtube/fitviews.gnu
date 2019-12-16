#!/bin/bash

set -eu

gnu=/tmp/$$.gnu

dat=${1:-data/views.log}

cat > $gnu << 'EOF'
set xdata time
set timefmt "%s"
set format x "%a\n%b%d\n%H:%M"

hour=3600

set mouse mouseformat 4
set key bottom right
x0=system("awk 'END {print $1}' ".dat)

set fit quiet
set fit logfile '/dev/null'

f(x) = fa*(x-x0) + fb
ft=24
fit [x0-hour*ft:x0] f(x) dat via fa,fb
fname=sprintf("%.2g views/hour (last %g hours)", fa*hour, ft)

g(x) = ga*(x-x0) + gb
gt=24*7
fit [x0-hour*gt:x0] g(x) dat via ga,gb
gname=sprintf("%.2g views/hour (last %g hours)", ga*hour, gt)

plot dat u 1:2 w l lw 2 ,\
  f(x) t fname ,\
  g(x) t gname
EOF

gnuplot -e "dat='$dat'" $gnu -

rm "$gnu"
