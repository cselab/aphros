SHELL := /bin/bash

CC = mpicxx
LD = mpicxx

SHELL := /bin/bash

# core
bs ?= 16
ap ?= float
omp ?= 1
align ?= 16

hdf ?= 1

CPPFLAGS+= $(extra)

ifeq "$(config)" "release"
	OPTFLAGS += -O3 -fno-expensive-optimizations -falign-functions=16 -DNDEBUG
	OPTFLAGS += -DNDEBUG
else
	CPPFLAGS += -g
endif

ifeq "$(omp)" "1"
	CPPFLAGS += -fopenmp
	OPTFLAGS += -fopenmp
endif

CPPFLAGS += -D_ALIGNBYTES_=$(align) -D_BLOCKSIZE_=$(bs) -D_BLOCKSIZEX_=$(bs) -D_BLOCKSIZEY_=$(bs) -D_BLOCKSIZEZ_=$(bs) 

ifeq "$(hdf)" "1"
	CPPFLAGS += -D_USE_HDF_
	LIBS += -lhdf5
endif

LIBS += -lstdc++ -lm -lz



CPPFLAGS += -I. -ICubism/

# hypre
include ~/bin/hypre.prefix
HYPRE_DIR = $(P)
LIBS += -L$(HYPRE_DIR)/lib -lHYPRE -lm
CPPFLAGS += -I$(HYPRE_DIR)/include

VPATH := Cubism/source/ CubismDistr hydro

.DEFAULT_GOAL := mfer

OBJECTS = main.o
OBJECTS += Cubism.o
OBJECTS += Local.o

all: mfer run

run:
	cd ../sim/sim01 && ./job

.PHONY: run

mfer: $(OBJECTS)
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o mfer

cleanall:
	make clean;

