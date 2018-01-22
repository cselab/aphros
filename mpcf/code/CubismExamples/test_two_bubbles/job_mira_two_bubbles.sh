#!/bin/bash
date; runjob --np 64 --ranks-per-node 1 --cwd $PWD --envs OMP_NUM_THREADS=64 --envs XLSMPOPTS=parthds=64 --envs PAMI_DEVICE=B --envs BG_MEMSIZE=16384 --envs BG_THREADLAYOUT=2 --envs OMP_STACKSIZE=4M --envs OMP_SCHEDULE="dynamic,1" --envs PAMID_COLLECTIVES=1 --envs PAMI_MEMORY_OPTIMIZED=1 --envs BG_SHAREDMEMSIZE=512 --envs BG_MAPCOMMONHEAP=0 --envs BG_SMP_FAST_WAKEUP=YES --envs L1P_POLICY="dcbt" --envs L1P_DEPTH=2 --envs PAMID_THREAD_MULTIPLE=1 --envs PAMID_VERBOSE=1 --envs PAMID_MAX_COMMTHREADS=1 --envs OMP_WAIT_POLICY=PASSIVE --envs OMP_PROC_BIND=FALSE --envs USEMAXTHREADS=0 --envs MYROUNDS=1 --envs DARSHAN_DISABLE=1 --block $COBALT_PARTNAME ${COBALT_CORNER:+--corner} $COBALT_CORNER ${COBALT_SHAPE:+--shape} $COBALT_SHAPE : ../../../CUBISM-MPCF/CubismApps/MPCFcluster/makefiles/mpcf-cluster \
-sim cloud -restart 0 \
-bpdx 4 -bpdy 4 -bpdz 4 -xpesize 4 -ypesize 4 -zpesize 4 -extent 1.0 \
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
