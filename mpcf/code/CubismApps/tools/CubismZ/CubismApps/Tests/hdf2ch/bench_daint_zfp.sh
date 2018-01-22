#2.09 3.039666e-08 4.834530e-12 3.086598e-10 5.172905e-10 2.881949e-12      15.3114     184.6956
#2.09 3.039666e-08 4.834530e-12 3.086598e-10 5.172905e-10 2.881949e-12      15.3114     184.6956
#2.14 3.039666e-08 4.034440e-11 2.575781e-09 1.453816e-09 8.099561e-12      14.9413     175.7201
#2.68 1.823799e-07 2.110281e-08 1.347305e-06 4.463491e-08 2.486718e-10      11.9413     145.9768
#4.03 2.431732e-06 5.146198e-07 3.285581e-05 6.613067e-07 3.684298e-09       7.9413     122.5622
#6.48 1.860275e-05 4.125325e-06 2.633807e-04 5.287230e-06 2.945643e-08       4.9408     104.5057
#16.74 1.378184e-04 2.225498e-05 1.420865e-03 2.807824e-05 1.564306e-07       1.9112      90.0029
#32.50 1.023577e-03 8.803176e-05 5.620373e-03 1.176064e-04 6.552133e-07       0.9847      77.5617

#make clean;make CC=CC wavz=1 zlib=1
make clean;make CC=CC zfp=1

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
