bs=32
ds=512
nb=$(echo "$ds/$bs" | bc)

make clean

# stage-0



# stage-1
make all wavz=1 zlib=1

# stage-2
rm -f output00000.StreamerGridPointIterative.channel0
./dat2ch  -bpdx $nb -bpdy $nb -bpdz $nb -sim io -threshold $1 -wtype_write 1

# stage-3
./ch2dat -simdata ./output00000.StreamerGridPointIterative.channel0 -datfile ch0 -wtype_read 1

# stage-4
./compdat ch0.dat ref.dat ./output00000.StreamerGridPointIterative.channel0




