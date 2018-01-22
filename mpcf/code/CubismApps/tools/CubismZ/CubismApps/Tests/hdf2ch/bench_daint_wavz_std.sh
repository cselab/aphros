#RES:         1.64 1.053548e-04 3.807098e-07 2.430635e-05 1.166185e-06 6.497098e-09     19.4798     117.6349
#RES:         2.28 2.009584e-03 1.541890e-05 9.844173e-04 3.327892e-05 1.854048e-07     14.0658      88.5269
#RES:        23.00 3.922308e-02 1.465144e-04 9.354185e-03 4.526027e-04 2.521558e-06	 1.3916      65.8559

#RES:         1.86 2.723540e-05 1.495959e-06 9.550924e-05 2.039599e-06 1.136309e-08     17.2091     112.7794
#RES:        12.30 3.553369e-04 2.263720e-05 1.445268e-03 2.920431e-05 1.627042e-07	 2.6026      89.6614
#RES:        60.30 2.880630e-03 8.630487e-05 5.510120e-03 1.240357e-04 6.910325e-07	 0.5307      77.0994


#make clean;make CC=CC wavz=1 zlib=1
make clean;make CC=CC wavz=1 zstd=1

#h5file=../../../../fabdata/data_010000-p.h5
#wt=1
#th=0.01

h5file=$1
wt=$2
th=$3

export OMP_PROC_BIND=TRUE;
export OMP_NUM_THREADS=1;  srun --ntasks=1 -c 12 --threads-per-core=1 ./hdf2ch -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $h5file  -outdata tmp  -threshold $th -wtype_write $wt
export OMP_NUM_THREADS=2;  srun --ntasks=1 -c 12 --threads-per-core=1 ./hdf2ch -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $h5file  -outdata tmp  -threshold $th -wtype_write $wt
export OMP_NUM_THREADS=4;  srun --ntasks=1 -c 12 --threads-per-core=1 ./hdf2ch -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $h5file  -outdata tmp  -threshold $th -wtype_write $wt
export OMP_NUM_THREADS=6;  srun --ntasks=1 -c 12 --threads-per-core=1 ./hdf2ch -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $h5file  -outdata tmp  -threshold $th -wtype_write $wt
export OMP_NUM_THREADS=8;  srun --ntasks=1 -c 12 --threads-per-core=1 ./hdf2ch -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $h5file  -outdata tmp  -threshold $th -wtype_write $wt
export OMP_NUM_THREADS=10; srun --ntasks=1 -c 12 --threads-per-core=1 ./hdf2ch -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $h5file  -outdata tmp  -threshold $th -wtype_write $wt
export OMP_NUM_THREADS=12; srun --ntasks=1 -c 12 --threads-per-core=1 ./hdf2ch -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $h5file  -outdata tmp  -threshold $th -wtype_write $wt
#export OMP_NUM_THREADS=16; srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $h5file  -outdata tmp  -threshold $th -wtype_write $wt
#export OMP_NUM_THREADS=20; srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $h5file  -outdata tmp  -threshold $th -wtype_write $wt
export OMP_NUM_THREADS=24; srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $h5file  -outdata tmp  -threshold $th -wtype_write $wt
