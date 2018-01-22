#!/bin/bash

#BSUB -J cbc             # provide a job name
#BSUB -o cbc_%J.out      # specify output file
#BSUB -e cbc_%J.err      # specify error file
#BSUB -W 00:30                     # wall clock time
#BSUB -n 1                         # number of processors
#BSUB -R "span[ptile=48]"          # number of processors per node

export OMP_NUM_THREADS=1;

# TODO: (fabianw; Tue 08 Sep 2015 11:23:12 AM CEST) WARNING: I have changed the
# parameter nu1 to mu1.  I do not know wether the set value corresponds to nu
# or mu, the solver expects mu -> please check the value before running this.

# mpirun -np 1 --npernode 1 --cpus-per-proc 48 \
# ../../CubismApps/MPCFcluster/makefiles/mpcf-cluster \
# -sim hit -restart 0 \
# -bpdx 1 -bpdy 1 -bpdz 1 -xpesize 1 -ypesize 1 -zpesize 1 -extent 6.28319 \
# -tend 2.5 -cfl 0.3 -nsteps 1 \
# -bcperiodicx 1 -bcperiodicy 1 -bcperiodicz 1 \
# -hllc 1 -mollfactor 2 -state 1 \
# -rho1 1.0 -p1 179.44 -mu1 0.0007431 -g1 1.4 -pc1 0.0 -g2 1.4 -pc2 0.0 \
# -io 1 -analysisperiod 10 -saveperiod 500 -dumpperiod 100 \
# -verbosity 1 \
# -kernels cpp \
# -dispatcher omp
