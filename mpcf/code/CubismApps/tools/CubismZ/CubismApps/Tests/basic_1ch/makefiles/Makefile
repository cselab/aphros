SHELL := /bin/bash

CC = mpicxx
LD = mpicxx
bs ?= 32
ap ?= float
config=release

include ../../../Makefile.config

mpi-inc ?=.
mpi-lib ?=.

ifneq "$(findstring silvio,$(shell hostname))" ""
	mpi-inc = /home/menahel/openmpi-1.5.5/include/
	mpi-lib = /home/menahel/openmpi-1.5.5/lib/
endif

ifneq "$(findstring brutus,$(shell hostname))" ""
	mpi-inc = ${MPI_ROOT}/include/
	mpi-lib = ${MPI_ROOT}/lib/
endif

ifneq "$(findstring rosa,$(shell hostname))" ""
	mpi-inc = $(CRAY_MPICH2_DIR)/include/
	mpi-lib = $(CRAY_MPICH2_DIR)/lib/
endif

ifneq "$(findstring daint,$(shell hostname))" ""
        mpi-inc = $(CRAY_MPICH2_DIR)/include/
        mpi-lib = $(CRAY_MPICH2_DIR)/lib/
endif

CPPFLAGS += -I../../../Cubism/source/ -I../../../Compressor/source/ -I../source
CPPFLAGS +=  -I$(mpi-inc)

#CPPFLAGS += -I${HOME}/usr/fpc 
#LIBS += -L${HOME}/usr/fpc -lfpc


ifneq "$(findstring rosa,$(shell hostname))" ""
	ifeq "$(CC)" "icc"
		LIBS += -L$(mpi-lib) -lmpich_intel -lmpichcxx_intel
	else
		LIBS += -L$(mpi-lib) -lmpich -lmpichcxx
	endif
else ifneq "$(findstring daint,$(shell hostname))" ""
        ifeq "$(CC)" "icc"
                LIBS += -L$(mpi-lib) -lmpich_intel -lmpichcxx_intel
        else
                LIBS += -L$(mpi-lib) -lmpich -lmpichcxx
        endif
else
	ifeq "$(findstring mpi,$(CC))" ""
		LIBS += -L$(mpi-lib) -lmpi -lmpi_cxx
	endif
endif

##################
VPATH := ../../../Cubism/source/ ../../../Compressor/source/ ../source/
.DEFAULT_GOAL := tests

OBJECTS = main.o WaveletCompressor.o

all: tests data_generator

tests: $(OBJECTS)
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

data_generator: data_generator.cpp
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o tests data_generator out* in.dat
