#!/bin/sh -u

cd $REPO/examples/205_multivof/bidisperse

cat > par.py << EOF
Lwide = 0.2
Lnarrow = 0.15
Hwide = 0.2
np = 8
nx = 64
tmax = 0.4
dumplist = ""
dumppoly = 0
dump_field_dt = 0.04
dump_traj_dt = 0.04
EOF

make run

(cd vis &&  ./vis.py --bubble_slice ../sm_*.vtk)

rsync -av "$(pwd)" /results/
