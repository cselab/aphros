#!/bin/sh

./age -p %.csv -f cl /u/fall/traj_*0.csv
./color -p '%.vtk' -f age -k cl *0.csv -- /u/fall/sm_*0.vtk
