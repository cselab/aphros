# Clean everything
make cleanall

# Compile
#make wavz=0 drain=1 zlib=0
make wavz=1 drain=0 zlib=1

# Run with 1 MPI process (only option for this case)
./tests -sim io -bpdx 4 -bpdy 2 -bpdz 2 

# Convert to HDF5
../../../tools/ch2hdf/ch2hdf -simdata output00000.StreamerGridPointIterative.channel0 -h5file out_ch0 
../../../tools/ch2hdf/ch2hdf -simdata output00000.StreamerGridPointIterative.channel2 -h5file out_ch2

