#!/bin/bash
set -x #echo on

##################
h5file=$1

./genref.sh $h5file

#./build.sh 11  # make all sz=1
#./bench_sz1.sh ./test_sz.sh $h5file

#./build.sh 8   # make all fpzip=1 
#./bench_fpzip1.sh ./test_fpzip.sh $h5file

#./build.sh 2   # make all wavz=1 lzma=1 shuffle3=1
#./bench_wavz1.sh ./test_wavz1.sh $h5file

#exit

./build.sh 2   # make all wavz=1 lzma=1 shuffle3=1
./bench_wavz.sh ./test_wavz.sh $h5file

./build.sh 14   # make all wavz=1 lzm=1 shuffle3=1
./bench_wavz.sh ./test_zstd.sh $h5file

exit

./build.sh 8   # make all fpzip=1 fpzip_eps=24
./bench_fpzip24.sh ./test_fpzip.sh $h5file

./build.sh 9   # make all fpzip=1 fpzip_eps=16
./bench_fpzip16.sh ./test_fpzip.sh $h5file

fpz24=`echo $h5file | sed 's/\// /g' | awk '{print $NF}'`_fpzip24_res.txt
fpz16=`echo $h5file | sed 's/\// /g' | awk '{print $NF}'`_fpzip16_res.txt
fpz=`echo $h5file | sed 's/\// /g' | awk '{print $NF}'`_fpzip_res.txt
cat $fpz24 $fpz16 > $fpz


./build.sh 10  # make all zfp=1
./bench_zfp.sh ./test_zfp.sh $h5file

./build.sh 11  # make all sz=1
./bench_sz.sh ./test_sz.sh $h5file
