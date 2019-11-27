.PHONY: all build run clean

all: build run

build: subcomm.cpp
	mpic++ -g -fopenmp -O2 subcomm.cpp -o main

run: build
	mpirun -n 2 ./main

clean:
	rm -v main
