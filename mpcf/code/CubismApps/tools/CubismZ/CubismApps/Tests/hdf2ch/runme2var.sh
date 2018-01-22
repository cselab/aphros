bs=32
ds=512
nb=$(echo "$ds/$bs" | bc)
#h5file=data-301-g.h5
h5file=data-301-p.h5

#make clean
#make all wavz=1 zlib=1

rm tmp00000.StreamerGridPointIterative.channel0

./hdf2ch -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $h5file  -outdata tmp  -threshold $1 -wtype_write 3

mpirun -n 8 ./ch2diff -simdata1 tmp00000.StreamerGridPointIterative.channel0  -simdata2 ref.channel0 -wtype 3

