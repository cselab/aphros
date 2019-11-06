#!/bin/sh

set -eu

t=/tmp/age.$$.vtk
trap 'rm $t; exit 1' 1 2 3 4 15
for i
do
    echo >&2 "$i"
    <"$i" ./wall | ./remove nn l c > $t &&
        mv $t "$i"
done
rm -f $t
