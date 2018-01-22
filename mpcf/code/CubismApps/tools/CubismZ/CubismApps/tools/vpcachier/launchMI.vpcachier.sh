#!/bin/sh

LOCARGS="--block $COBALT_PARTNAME ${COBALT_CORNER:+--corner} $COBALT_CORNER ${COBALT_SHAPE:+--shape} $COBALT_SHAPE"
echo "Cobalt location args: $LOCARGS" >&2

#MIN=-54.67
#MAX=103.49

   runjob $LOCARGS \
  --verbose=INFO \
  --label \
  --np 64 \
  --ranks-per-node 1 \
  --cwd $PWD \
  --exe $PWD/vpcachier \
  --envs PAMI_DEVICE=B \
  --envs BG_MEMSIZE=16384 \
  --envs BG_THREADLAYOUT=2 \
  --envs XLSMPOPTS=parthds=64 \
  --envs PAMID_COLLECTIVES=1 \
  --envs PAMI_MEMORY_OPTIMIZED=1 \
  --envs BG_SHAREDMEMSIZE=512 \
  --envs BG_MAPCOMMONHEAP=0 \
  --args -simdata datawavelet000050.StreamerGridPointIterative.channel5 \
  --args -vp vpcache1.little \
  --args -swapW

#  --args -vp vpcache1.big \

#  --envs MINVAL=$MIN \
#  --envs MAXVAL=$MAX

#  --args -min $MIN  \
#  --args -max $MAX

#  --args -restart 0

#  --mapping TABCDE \
#  --envs BG_COREDUMPONEXIT=1 \

#  --envs BGLOCKLESSMPIO_F_TYPE=0x47504653


#  --envs PAMID_CONTEXT_POST=
#  --envs MUSPI_RECFIFOSIZE=67108864 \
#  --envs PAMID_CONTEXT_POST=1 \
