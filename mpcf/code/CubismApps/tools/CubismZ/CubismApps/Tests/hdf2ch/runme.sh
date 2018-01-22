bs=32
ds=512
nb=$(echo "$ds/$bs" | bc)
#h5file=data-301-g.h5
h5file=data-301-p.h5
#h5file=ch4s300.h5
#h5file=ch5s300.h5

make clean

#make all wavz=1 zlib=1
#make all wavz=1 lzma=1
#make all wavz=1 lzma=1 shuffle3=1
#make all wavz=1 lzma=1 shuffle3=1 zerobits=4
#make all wavz=1 lzma=1 shuffle3=1 zerobits=8
#make all wavz=1 lzma=1 shuffle3=1 zerobits=12
#make all wavz=1 lzma=1 shuffle3=1 zerobits=16
#make all fpzip=1 fpzip_eps=32
#make all fpzip=1 fpzip_eps=24
#make all fpzip=1 fpzip_eps=16
make all zfp=1

#make all sz=1
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/jvm/java-1.8.0/jre/lib/amd64/server

#make all isa=1


rm tmp00000.StreamerGridPointIterative.channel0

./hdf2ch -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $h5file  -outdata tmp  -threshold 0.008 -wtype_write 3

mpirun -n 8 ./ch2diff -simdata1 tmp00000.StreamerGridPointIterative.channel0  -simdata2 ref.channel0 -wtype 3

