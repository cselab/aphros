#!/bin/bash
set -x #echo on

##################
h5file=$1

./genref.sh $h5file

./build.sh 8   # make all fpzip=1 fpzip_eps=24
./bench_fpzip24.sh ./test_fpzip.sh $h5file

./build.sh 9   # make all fpzip=1 fpzip_eps=16
./bench_fpzip16.sh ./test_fpzip.sh $h5file

fpz24=`echo $h5file | sed 's/\// /g' | awk '{print $NF}'`_fpzip24_res.txt
fpz16=`echo $h5file | sed 's/\// /g' | awk '{print $NF}'`_fpzip16_res.txt
fpz=`echo $h5file | sed 's/\// /g' | awk '{print $NF}'`_fpzip_res.txt
cat $fpz24 $fpz16 > $fpz
