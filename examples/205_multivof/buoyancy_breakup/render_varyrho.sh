#!/bin/sh -eu

cd vis
./single.py \
  --files0 ../rho0.00125/sm_*.vtk \
  --files1 ../rho0.01000/sm_*.vtk \
  "$@"
