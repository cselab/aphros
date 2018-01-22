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
#h5file=data-301-g.h5
#h5file=data-301-p.h5
#h5file=ch4s300.h5
#h5file=ch5s300.h5

nb=2
rm -f ref.channel0

make clean
make all
./hdf2ch -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $h5file  -outdata c1 
mv c100000.StreamerGridPointIterative.channel0 ref.channel0




