export OMP_PROC_BIND=TRUE;
export OMP_NUM_THREADS=1;  srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=2;  srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=4;  srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=6;  srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=8;  srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=10; srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=12; srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=14; srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=16; srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=18; srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=20; srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=22; srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
export OMP_NUM_THREADS=24; srun --ntasks=1 -c 12 --threads-per-core=2 ./hdf2ch -xpesize 1 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata ../../../../fabdata/data_010000-p.h5 -threshold 0.0001 -outdata c1 
