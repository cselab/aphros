#!/bin/bash

# some hints
# - use at least 3 blocks per direction
# - specify -kernels cpp
# - provide path to 'mpcf-cluster'
./mpcf-cluster -sim cloud -restart 0 \
-bpdx 3 -bpdy 3 -bpdz 3 -xpesize 1 -ypesize 1 -zpesize 1 -extent 0.6 \
-tend 10.0 -cfl 0.3 -nsteps 0 \
-sponge 1 \
-hllc 1 -mollfactor 2 -state 1 \
-pref 100.0 \
-g1 6.59 -pc1 4.049e3 -g2 1.4 -pc2 1.0 \
-io 1 -analysisperiod 10 -saveperiod 100 \
-vp 1 -dumpdt 0.02 -vpeps 1e-3 -vpchannels 0123456 \
-verbosity 1 \
-kernels cpp \
-dispatcher omp 
