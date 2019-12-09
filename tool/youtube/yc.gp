#!/bin/sh

awk '
{
    month = $2
    day = $3
    split($4, tmp, ":")
    hour = tmp[1]
    minut = tmp[2]
    year = $6
    count = $7
    label = day "/" hour ":" minut
    epoch = (minut + 60*(hour + 24*day))
    print epoch, count, label
}' "$@" > /tmp/t

gnuplot -p <<'!'
set yrange [0:1000]
set key off
plot "/tmp/t" u 1:2 w l, \
     "<awk '!(NR % 50)' /tmp/t"       u 1:2:xticlabels(3) w p
!
