export OMP_PROC_BIND=TRUE;
#export OMP_NUM_THREADS=1; srun --ntasks=1 -c 1  ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=2; srun --ntasks=1 -c 2  ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=4; srun --ntasks=1 -c 4  ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=6; srun --ntasks=1 -c 6  ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=8; srun --ntasks=1 -c 8  ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=10; srun --ntasks=1 -c 10  ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=12; srun --ntasks=1 -c 12  ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
