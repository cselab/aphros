#!/bin/bash

#SBATCH --job-name=simple_cloud        # job name
#SBATCH --ntasks=64                    # number of tasks
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00                # rum time hours:minutes:seconds
#SBATCH --output=simple_cloud_jobID-%j.out    # output file name
#SBATCH --error=simple_cloud_1_jobID-%j.err   # error file name
#SBATCH --account=s500                 # project id (IMPORTANT)

export OMP_NUM_THREADS=8;

aprun -n $SLURM_NTASKS -N $SLURM_NTASKS_PER_NODE -d $OMP_NUM_THREADS \
../../../CUBISM-MPCF/CubismApps/MPCFcluster/makefiles/mpcf-cluster \
-sim cloud -restart 0 \
-bpdx 6 -bpdy 6 -bpdz 6 -xpesize 4 -ypesize 4 -zpesize 4 -extent 1.0 \
-tend 10.0 -cfl 0.3 -nsteps 0 \
-sponge 1 \
-hllc 1 -mollfactor 2 -state 1 \
-pref 100.0 \
-g1 6.59 -pc1 4.049e3 -g2 1.4 -pc2 1.0 \
-io 1 -analysisperiod 10 -saveperiod 100 \
-vp 1 -dumpdt 0.02 -vpeps 1e-3 -vpchannels 0123456wcMm \
-verbosity 1 \
-kernels cpp \
-dispatcher omp -ncores 8
