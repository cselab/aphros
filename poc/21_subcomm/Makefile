.PHONY: all build run clean

CC = mpic++

all: build run

build: subcomm.cpp
	$(CC) -g -fopenmp -O2 subcomm.cpp -o main

run: build
	mpirun -n 2 ./main

daint: build
	export MPICH_ENV_DISPLAY=1
	srun -n 24 -N 2 --ntasks-per-node=12 --cpus-per-task=1 ./main

clean:
	rm -v main
