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

CPPFLAGS += -I../../../Cubism/source/ -I../../../Compressor/source/ -I. 
CPPFLAGS += -I../../../../CubismApps/tools/reader/
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
VPATH := ../../../Cubism/source/ ../../../Compressor/source/ .
.DEFAULT_GOAL := dat2ch 

OBJECTS1 = dat2ch.o WaveletCompressor.o
OBJECTS2 = ch2dat.o WaveletCompressor.o

all: dat2ch data_generator ch2dat compdat

#ch2dat: ch2dat.o WaveletCompressor.o
#	$(CC) $(CPPFLAGS) WaveletCompressor.o  ch2dat.o  -o ch2dat $(LIBS) -ldl

#ch2dat.o: ch2dat.cpp
#       $(CC) $(CPPFLAGS) -c ch2dat.cpp -o ch2dat.o

dat2ch: $(OBJECTS1)
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

ch2dat: $(OBJECTS2)
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

data_generator: data_generator.cpp
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

compdat: compdat.c
	$(CC) $(OPTFLAGS) $(extra) compdat.c -o compdat


%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o dat2ch data_generator ch2dat compdat out* 

cleanall:
	rm -f *.o dat2ch data_generator ch2dat compdat out* in.dat


