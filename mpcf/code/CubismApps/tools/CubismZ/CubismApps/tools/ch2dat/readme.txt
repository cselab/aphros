# module load gcc/4.7 open_mpi/1.6.2 hdf5_para
# make
# mpirun -mca btl sm,self -np 4 vp2hdf -simdata ../data/datawavelet00000.StreamerGridPointIterative.channel4 -h5file h5outputfile


