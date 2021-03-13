#!/bin/sh -u

cd $REPO/examples/205_multivof/waterfall

cat > par.py << EOF
np = 8
nn = 32
tmax = 2
dumpdt = 0.2
EOF

cat > add.conf << EOF
set string dumplist
set int dumppoly 0
EOF

make run

(cd vis && ./viswhite.py --draft --ray 0 --ambient 0 ../sm_*.vtk)

rsync -av "$(pwd)" /results/
