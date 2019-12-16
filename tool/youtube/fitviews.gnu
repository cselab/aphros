#!/bin/bash

gnu=/tmp/$$.gnu
dat=/tmp/$$.dat

rsync "falcon:views.log" "$dat"

cat > $gnu << 'EOF'
set xdata time
set timefmt "%s"
set format x "%a\n%b%d\n%H:%M"
hour=3600
set xtics hour*24
set mouse mouseformat 4
file="views.log"
x0=system("date +%s")

set fit quiet
set key bottom right

f(x) = fa*(x-x0) + fb
fit [x0-hour*24:x0] f(x) dat via fa,fb
fname=sprintf("%.2g views/hour", fa*hour)

g(x) = ga*(x-x0) + gb
fit [x0-hour*24*7:x0] g(x) dat via ga,gb
gname=sprintf("%.2g views/hour", ga*hour)

plot dat u 1:2 w l lw 2 ,\
  f(x) t fname ,\
  g(x) t gname
EOF

gnuplot -e "dat='$dat'" $gnu -

rm "$dat" "$gnu"
