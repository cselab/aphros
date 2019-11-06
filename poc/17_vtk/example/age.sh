#!/bin/sh

set -eu

./age -p %.csv -f cl /u/fall/traj_*.csv
./color -p a.%.vtk -f age -k cl *0.csv -- /u/fall/sm_*0.vtk

