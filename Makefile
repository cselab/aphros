SHELL := /bin/bash

CC = mpic++
LD = mpic++

bsx ?= 16
bsy ?= $(bsx)
bsz ?= 16
ap ?= float
omp ?= 1
align ?= 16
hdf ?= 1

CPPFLAGS = -I./Cubism/source/
CPPFLAGS += -D_ALIGNBYTES_=$(align) -D_BLOCKSIZEX_=$(bsx) -D_BLOCKSIZEY_=$(bsy) -D_BLOCKSIZEZ_=$(bsz)
LIBS = -lstdc++ -lm -lz

ifeq "$(ap)" "float"
	CPPFLAGS += -D_FLOAT_PRECISION_
endif

ifeq "$(omp)" "1"
	CPPFLAGS += -fopenmp
endif


ifeq "$(hdf)" "1"
	CPPFLAGS += -I/opt/hdf5_openmpi/include -D_USE_HDF_
	LIBS += -L/opt/hdf5_openmpi/lib -lhdf5
endif

LIBS += -ldl

.DEFAULT_GOAL := main

all: test-cubismo

main: main.cpp
	$(CC) $(CPPFLAGS) $(extra) $^ -o $@ $(LIBS)

clean:
	rm -f *.o main 

cleandeep: clean
	find . -iname "*~" -exec rm -f {} \;

run:
	mpirun -np 2 ./main
