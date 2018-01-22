#!/bin/sh
NODES=`cat $COBALT_NODEFILE | wc -l`
PROCS=$((NODES * 16))

mpirun \
-f $COBALT_NODEFILE \
-n 1 \
/home/chatzidp/mywork/CUBISM-MPCF/CubismApps/tools/vpcachier/vpcachier \
-simdata datawavelet00050.StreamerGridPointIterative.channel5 \
-vp vpcache1_tk \
-swap
