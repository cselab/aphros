SHELL := /bin/bash

CC = mpic++
bgq ?= 0
omp ?= 1
ap ?= float
fftw ?= 1
hdf ?= 1
config ?= release

ifeq "$(omp)" "1"
	ifeq "$(CC)" "icc"
		CPPFLAGS += -openmp
		OPTFLAGS += -openmp
	else
		ifeq "$(bgq)" "1"
			CPPFLAGS += -qsmp=speculative
			OPTFLAGS += -qsmp=speculative
		else
			CPPFLAGS += -fopenmp
			OPTFLAGS += -fopenmp
		endif
	endif
endif

ifeq "$(ap)" "float"
	CPPFLAGS += -D_FLOAT_PRECISION_
endif

ifeq "$(config)" "release"
	ifeq "$(CC)" "icc"
		OPTFLAGS+= -DNDEBUG -O3 -xHOST  -ip -ansi-alias -fno-fnalias -inline-level=1
		ifeq "$(ap)" "float"
			CPPFLAGS += -fp-model precise
		else
			CPPFLAGS += -fast
		endif
	else
		ifeq "$(bgq)" "1"
			OPTFLAGS += -O3 -qsmp=speculative -qhot=level=1  -qarch=qp -qtune=qp -qprefetch=aggressive
		else
			OPTFLAGS += -O3 -fno-expensive-optimizations -DNDEBUG
		endif
	endif
endif

hdf-inc ?=.
hdf-lib ?=.
fftw-inc ?=.
fftw-lib ?=.

ifneq "$(findstring brutus,$(shell hostname))" ""
	mpi-inc = ${MPI_ROOT}/include/
	mpi-lib = ${MPI_ROOT}/lib/
        hdf-inc = $(HDF5_ROOT)/include/
        hdf-lib = $(HDF5_ROOT)/lib/
endif

ifneq "$(findstring fen,$(shell hostname))" ""
       hdf-inc = ${HDF5_INC}
       hdf-lib = ${HDF5_LIB}
       fftw-inc = ${FFTW3_INC}
       fftw-lib = ${FFTW3_LIB}
endif

ifneq "$(findstring cetus,$(shell hostname))" ""
       hdf-inc = /soft/libraries/hdf5/current/cnk-xl/current/include
       hdf-lib = /soft/libraries/hdf5/current/cnk-xl/current/lib
       fftw-inc = /soft/libraries/alcf/current/xl/FFTW3/include/
       fftw-lib = /soft/libraries/alcf/current/xl/FFTW3/lib/
endif

ifneq "$(findstring mira,$(shell hostname))" ""
       hdf-inc = /soft/libraries/hdf5/current/cnk-xl/current/include
       hdf-lib = /soft/libraries/hdf5/current/cnk-xl/current/lib
       fftw-inc = /soft/libraries/alcf/current/xl/FFTW3/include/
       fftw-lib = /soft/libraries/alcf/current/xl/FFTW3/lib/
endif

ifneq "$(findstring juq,$(shell hostname))" ""
         hdf-inc = ${HDF5_INCLUDE}
         hdf-lib = ${HDF5_LIB}
         fftw-inc = ${FFTW_INCLUDE}
         fftw-lib = ${FFTW_LIB}
endif

ifeq "$(fftw)"  "1"
#FFTW always in double precision
	CPPFLAGS += -I$(fftw-inc) -D_USE_FFTW_
#	LIBS += -L$(fftw-lib) -lfftw3 -lfftw3_omp -lfftw3_mpi
        ifeq "$(bgq)" "1"
              LIBS += -L$(fftw-lib) -lfftw3_omp -lfftw3_mpi -lfftw3
        else
              LIBS += -L$(fftw-lib) -lfftw3 -lfftw3_omp -lfftw3_mpi
        endif
endif

ifeq "$(hdf)" "1"
	CPPFLAGS += -I$(hdf-inc) -D_USE_HDF_
	LIBS += -L$(hdf-lib) -lhdf5

	ifneq "$(findstring juq,$(shell hostname))" ""
	      LIBS += -L/bgsys/local/szip/lib/ -lsz -L/bgsys/local/zlib/lib/ -lz
	endif

        ifneq "$(findstring fen,$(shell hostname))" ""
              #CPPFLAGS += -DH5_USE_16_API
              #LIBS += -L$SZIP_LIB -lsz -L$ZLIB_LIB -lz
              LIBS += -L/cineca/prod/libraries/szip/2.1/bgq-xl--1.0/lib -lsz -L/cineca/prod/libraries/zlib/1.2.7/bgq-gnu--4.4.6/lib -lz
              CPPFLAGS += -D_ON_FERMI_
        endif

        ifneq "$(findstring cetus,$(shell hostname))" ""
              LIBS += -L/soft/libraries/alcf/current/xl/SZIP/lib/ -L/soft/libraries/alcf/current/xl/ZLIB/lib/ -lz
        endif

        ifneq "$(findstring mira,$(shell hostname))" ""
              LIBS += -L/soft/libraries/alcf/current/xl/SZIP/lib/ -L/soft/libraries/alcf/current/xl/ZLIB/lib/ -lz
        endif
endif

CPPFLAGS += -g

CPPFLAGS += -I../../../../Cubism/source/
# CPPFLAGS += -I../../Cubism/source -I../Cubism/source -I../../../Cubism/source/ -I../../MPCFcluster/source/ -I../../MPCFnode/source -I../../MPCcore/source -I.  -I../reader

##################
# VPATH := ../../../Cubism/source/ ../../MPCFcluster/source/ ../../MPCFnode/source/ .
.DEFAULT_GOAL := spectrum

OBJECTS = main.o ../FFTWBaseMPI.o eval_hit/evalHit.o

all: spectrum

spectrum: $(OBJECTS)
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o spectrum *~ ../*.o eval_hit/evalHit.o
