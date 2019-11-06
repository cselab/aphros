#!/bin/sh

set -eu

./age -p %.csv -f cl /u/fall/traj_*.csv
./color -p a.%.vtk -f age -k cl *0.csv -- /u/fall/sm_*0.vtk

t=/tmp/age.$$.vtk
trap 'rm $t; exit 1' 1 2 3 4 15
for i in a.*.vtk
do
    echo >&2 "$i"
    <$i | ./wall | ./remove nn l c > $t &&
        mv $t $i
done
rm -f $t
