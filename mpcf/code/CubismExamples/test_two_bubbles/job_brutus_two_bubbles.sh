#!/bin/bash

#BSUB -J two_bubbles             # provide a job name
#BSUB -o two_bubbles_%J.out      # specify output file
#BSUB -e two_bubbles_%J.err      # specify error file
#BSUB -W 08:00                   # wall clock time
#BSUB -n 384                      # number of processors
#BSUB -R "span[ptile=48]"        # number of processors per node

export OMP_NUM_THREADS=48;

mpirun -np 8 --npernode 1 --cpus-per-proc 48 \
../../../CUBISM-MPCF/CubismApps/MPCFcluster/makefiles/mpcf-cluster \
-sim cloud -restart 0 \
-bpdx 8 -bpdy 8 -bpdz 8 -xpesize 2 -ypesize 2 -zpesize 2 -extent 1.0 \
-tend 10.0 -cfl 0.3 -nsteps 0 \
-sponge 1 \
-hllc 1 -mollfactor 2 -state 1 \
-pref 100.0 \
-g1 6.59 -pc1 4.049e3 -g2 1.4 -pc2 1.0 \
-io 1 -analysisperiod 10 -saveperiod 100 \
-vp 1 -dumpdt 0.02 -vpeps 1e-3 -vpchannels 0123456wcMm \
-verbosity 1 \
-kernels cpp \
-dispatcher omp
