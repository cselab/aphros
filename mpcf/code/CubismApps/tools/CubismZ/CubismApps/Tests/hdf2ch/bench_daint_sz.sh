#for PARAM in 0.000001 0.000010 0.000100 0.001000 0.010000 0.100000 1 10 50 100 200 500 1000

#1.58 7.599164e-09 3.338764e-16 2.131628e-14 2.209016e-12 1.230696e-14      20.1981     232.0863
#1.75 7.599164e-08 1.153479e-08 7.364366e-07 3.483763e-08 1.940888e-10      18.2414     148.1293
#2.13 7.903130e-07 1.945383e-07 1.242026e-05 3.259476e-07 1.815932e-09      15.0523     128.7073
#2.88 7.963924e-06 3.445974e-06 2.200076e-04 4.923744e-06 2.743136e-08      11.0999     105.1244
#5.81 7.967723e-05 5.551320e-05 3.544231e-03 7.086710e-05 3.948176e-07       5.5039      81.9614
#10.87 7.968255e-04 6.125161e-04 3.910599e-02 7.583256e-04 4.224814e-06       2.9431      61.3732


#make clean;make CC=CC wavz=1 zlib=1
make clean;make CC=CC sz=1

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
#export OMP_NUM_THREADS=24; srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $h5file  -outdata tmp  -threshold $th -wtype_write $wt
