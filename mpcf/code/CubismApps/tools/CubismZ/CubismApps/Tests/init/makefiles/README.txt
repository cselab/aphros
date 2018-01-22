# compile for a specific configuration
make clean;make wavz=1 zlib=1 
#make clean;make wavz=0 drain=1 zlib=0 

# run this test
./tests -sim io -bpdx 4 -bpdx 4 -bpdz 4

# convert to HDF5 (ch2hdf must have been compiled as in the first step above) 
../../../tools/ch2hdf/ch2hdf -simdata output00000.StreamerGridPointIterative.channel0 -h5file out1

# visualize with Paraview 
paraview out1.xmf

