#!/bin/bash
set -x #echo on

##################
h5file=$1

./genref.sh $h5file

./build.sh 1   # make all wavz=1 zlib=1

./bench_wavz_type.sh ./test_wavz_type.sh $h5file 1
./bench_wavz_type.sh ./test_wavz_type.sh $h5file 2
./bench_wavz_type.sh ./test_wavz_type.sh $h5file 3
