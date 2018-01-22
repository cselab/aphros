#!/bin/bash
# @ job_type = bluegene
# @ job_name = simple_cloud
# @ output = $(job_name)_$(jobid).out
# @ error = $(job_name)_$(jobid).out
# @ environment = COPY_ALL
# @ wall_clock_limit = 00:10:00
# @ notification = never
# @ bg_size = 64
# @ bg_connectivity = MESH
# @ account_no = Pra09_2376
# @ queue

module load scalapack
runjob --np 64 --ranks-per-node 1 --envs OMP_NUM_THREADS=16 : ../../../CUBISM-MPCF/CubismApps/MPCFcluster/makefiles/mpcf-cluster \
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
-kernels qpx -pp 12.8 -pb 28.0 \
-dispatcher omp -ncores 16 -gsync 128 
