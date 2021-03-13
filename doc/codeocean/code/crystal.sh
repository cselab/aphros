#!/bin/sh -u

cd $REPO/examples/205_multivof/crystal

cat > par.py << EOF
np = 8
ny = 24
mu_liquid = 1e-3 * 0.25
sigma = 72e-3 * 0.05
rho_ratio = 0.1
erasevf_per_rel = 0.5
length = 3e-3
inletpressure_add = 0.5
tmax = 5.
lbuf = 0.2
erasevf_length = 0.15
dumplist = ""
dump_field_dt = 0.5
EOF

cat > add.conf << EOF
set int dumppoly 0
EOF

make run

(cd vis && ./vis.py --bubble_slice ../sm_*.vtk)

rsync -av "$(pwd)" /results/
