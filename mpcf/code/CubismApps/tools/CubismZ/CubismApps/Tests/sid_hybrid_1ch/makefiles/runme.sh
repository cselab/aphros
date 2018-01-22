# Clean everything
make cleanall

# Compile
#make wavz=0 drain=1 zlib=0 
make wavz=1 zlib=1

# Run with 1 MPI process
./tests -sim io -bpdx 4 -bpdy 2 -bpdz 2 
mv output00000.StreamerGridPointIterative.channel0 output00000.StreamerGridPointIterative.channel0_p1

# Run with 2 MPI processes
mpirun -n 2 ./tests -sim io -xpesize 2 -bpdx 2 -bpdy 2 -bpdz 2 
mv output00000.StreamerGridPointIterative.channel0 output00000.StreamerGridPointIterative.channel0_p2 

# Convert to HDF5
../../../tools/ch2hdf/ch2hdf -simdata output00000.StreamerGridPointIterative.channel0_p1 -h5file out_p1 
../../../tools/ch2hdf/ch2hdf -simdata output00000.StreamerGridPointIterative.channel0_p2 -h5file out_p2

# Compare
h5diff out_p1.h5 out_p2.h5 

