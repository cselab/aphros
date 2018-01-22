#!/bin/bash

export OMP_NUM_THREADS=48

JOBNAME='single_bubble_jc06'

bsub -n 384 -W 01:00 -R "span[ptile=48]" -J $JOBNAME -o ${JOBNAME}_%J.out -e ${JOBNAME}_%J.err \
    mpirun -np 8 --npernode 1 --cpus-per-proc 48 \
    ./mpcf-cluster -simConf single_bubble_jc06.conf "$@"
