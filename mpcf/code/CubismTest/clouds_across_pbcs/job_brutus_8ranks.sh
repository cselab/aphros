#!/bin/bash

#BSUB -J cbc             # provide a job name
#BSUB -o cbc_%J.out      # specify output file
#BSUB -e cbc_%J.err      # specify error file
#BSUB -W 00:20                     # wall clock time
#BSUB -n 384                        # number of processors
#BSUB -R "span[ptile=48]"          # number of processors per node

export OMP_NUM_THREADS=1;

# TODO: (fabianw; Tue 08 Sep 2015 11:24:23 AM CEST) WARNING: I have changed the
# parameter nu to mu.  I do not know wether the set value corresponds to nu
# or mu, the solver expects mu -> please check the value before running this.

# mpirun -np 8 --npernode 1 --cpus-per-proc 48 \
# /cluster/home04/mavt/ursular/2015_06_cubism_git/CUBISM-MPCF/CubismApps/MPCFcluster/makefiles/mpcf-cluster \
# -sim hitcloud -restart 0 \
# -bpdx 1 -bpdy 8 -bpdz 8 -xpesize 8 -ypesize 1 -zpesize 1 -extent 6.28319 \
# -tend 2.5 -cfl 0.3 -nsteps 1 \
# -bcperiodicx 1 -bcperiodicy 1 -bcperiodicz 1 \
# -hllc 1 -mollfactor 2 -state 1 \
# -rho1 1.0 -p1 1.396e4 -mu1 7.431e-4 -g1 6.59 -pc1 5.652e5 \
# -rho2 1.0e-3 -p2 3.267 -mu2 1.358e-5 -g2 1.4 -pc2 139.60 \
# -io 1 -analysisperiod 10000 -saveperiod 1000 -dumpperiod 10000 \
# -verbosity 1 \
# -kernels qpx \
# -dispatcher omp
