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


if [ -z "$3" ]
then
	wt=3	
else
	wt=$3
fi

#h5file=data-301-g.h5
#h5file=data-301-p.h5
#h5file=ch4s300.h5
#h5file=ch5s300.h5

#make clean
#make all wavz=1 zlib=1

rm tmp00000.StreamerGridPointIterative.channel0

./hdf2ch -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $h5file -outdata tmp  -threshold $2 -wtype_write $wt

mpirun -n 8 ./ch2diff -simdata1 tmp00000.StreamerGridPointIterative.channel0  -simdata2 ref.channel0 -wtype $wt

