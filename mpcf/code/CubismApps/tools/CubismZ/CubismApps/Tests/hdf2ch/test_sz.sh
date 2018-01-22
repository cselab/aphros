#!/bin/bash
set -x #echo on

if [ -z "$1" ]
then
    	echo "missing file"
        exit
else
    	h5file=$1
fi

bs=32
ds=512
nb=$(echo "$ds/$bs" | bc)
#wt=3
#h5file=data-301-g.h5
#h5file=data-301-p.h5
#h5file=ch4s300.h5
#h5file=ch5s300.h5

#make clean
#make all sz=1

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/jvm/java-1.8.0/jre/lib/amd64/server

rm tmp00000.StreamerGridPointIterative.channel0

./hdf2ch -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $h5file  -outdata tmp  -threshold $2

mpirun -n 8 ./ch2diff -simdata1 tmp00000.StreamerGridPointIterative.channel0  -simdata2 ref.channel0

