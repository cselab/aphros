#!/bin/sh -u

cd $REPO/examples/205_multivof/clustering

cat > add.conf << EOF
set string bubgen_path inline 0.015 0.007 0.015 0.002
set string list_path inline 0 0.03 0 1e5 0.015 1e5
set int dumppoly 0
set double dump_field_dt 0.01
set double tmax 0.1
set double bubgen_per 0.02
set double cflst 1
EOF

make m="32 24 32" bs="8 8 8" np=8 run

(cd vis && ./viswhite.py --draft --ray 0 --ambient 0 --camera exp ../sm_*.vtk)

rsync -av "$(pwd)" /results/
