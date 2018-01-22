#!/bin/bash
set -x #echo on

##################
h5file=/home/chatzidp/gitlab/fabdata/data_010000-p.h5

./genref.sh $h5file

./build.sh 2   # make all wavz=1 lzma=1 shuffle3=1
./bench_wavz.sh ./test_wavz.sh $h5file

./build.sh 8   # make all fpzip=1 fpzip_eps=24
./bench_fpzip24.sh ./test_fpzip.sh $h5file

./build.sh 9   # make all fpzip=1 fpzip_eps=16
./bench_fpzip24.sh ./test_fpzip.sh $h5file

fpz16=`echo $h5file | sed 's/\// /g' | awk '{print $NF}'`_fpzip16_res.txt
fpz24=`echo $h5file | sed 's/\// /g' | awk '{print $NF}'`_fpzip24_res.txt
fpz=`echo $h5file | sed 's/\// /g' | awk '{print $NF}'`_fpzip_res.txt
cat $fpz16 $fpz24 > $fpz


./build.sh 10  # make all zfp=1
./bench_zfp.sh ./test_zfp.sh $h5file

./build.sh 11  # make all sz=1
./bench_sz.sh ./test_sz.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/data_010000-rho.h5

./genref.sh $h5file

./build.sh 2   # make all wavz=1 lzma=1 shuffle3=1
./bench_wavz.sh ./test_wavz.sh $h5file

./build.sh 8   # make all fpzip=1 fpzip_eps=24
./bench_fpzip24.sh ./test_fpzip.sh $h5file

./build.sh 9   # make all fpzip=1 fpzip_eps=16
./bench_fpzip24.sh ./test_fpzip.sh $h5file

fpz16=`echo $h5file | sed 's/\// /g' | awk '{print $NF}'`_fpzip16_res.txt
fpz24=`echo $h5file | sed 's/\// /g' | awk '{print $NF}'`_fpzip24_res.txt
fpz=`echo $h5file | sed 's/\// /g' | awk '{print $NF}'`_fpzip_res.txt
cat $fpz16 $fpz24 > $fpz

./build.sh 10  # make all zfp=1
./bench_zfp.sh ./test_zfp.sh $h5file

./build.sh 11  # make all sz=1
./bench_sz.sh ./test_sz.sh $h5file


##################
h5file=/home/chatzidp/gitlab/fabdata/data_010000-a2.h5

./genref.sh $h5file

./build.sh 2   # make all wavz=1 lzma=1 shuffle3=1
./bench_wavz.sh ./test_wavz.sh $h5file

./build.sh 8   # make all fpzip=1 fpzip_eps=24
./bench_fpzip24.sh ./test_fpzip.sh $h5file

./build.sh 9   # make all fpzip=1 fpzip_eps=16
./bench_fpzip24.sh ./test_fpzip.sh $h5file

fpz16=`echo $h5file | sed 's/\// /g' | awk '{print $NF}'`_fpzip16_res.txt
fpz24=`echo $h5file | sed 's/\// /g' | awk '{print $NF}'`_fpzip24_res.txt
fpz=`echo $h5file | sed 's/\// /g' | awk '{print $NF}'`_fpzip_res.txt
cat $fpz16 $fpz24 > $fpz

./build.sh 10  # make all zfp=1
./bench_zfp.sh ./test_zfp.sh $h5file

./build.sh 11  # make all sz=1
./bench_sz.sh ./test_sz.sh $h5file
