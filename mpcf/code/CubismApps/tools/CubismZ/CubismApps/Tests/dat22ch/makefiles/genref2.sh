bs=32
ds=512
nb=$(echo "$ds/$bs" | bc)

make clean

# stage-0



# stage-1
make clean;make all	# no options, dummy compressions

# stage-2
./dat2ch  -bpdx $nb -bpdy $nb -bpdz $nb -sim io

# stage-3
./ch2dat -simdata ./output00000.StreamerGridPointIterative.channel0 -datfile ref

#  stage-4
#./compdat ch0.dat ref.dat ./output00000.StreamerGridPointIterative.channel0




