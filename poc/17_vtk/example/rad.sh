#!/bin/sh

set -eu

./rad -p r.%.vtk -f vf -k cl /u/fall/traj_*0.csv -- /u/fall/sm_*0.vtk
