SHELL := /bin/bash

CC ?= mpicxx
LD ?= mpicxx

include ../../Makefile.config

mpi-inc ?=.
mpi-lib ?=.
fftw-inc ?=.
fftw-lib ?=.

ifneq "$(findstring silvio,$(shell hostname))" ""
	mpi-inc = /home/menahel/openmpi-1.5.5/include/
	mpi-lib = /home/menahel/openmpi-1.5.5/lib/
	ifeq "$(CC)" "icc"
		fftw-inc = /home/hbabak/fftw-2.1.5_icc/include/
		fftw-lib = /home/hbabak/fftw-2.1.5_icc/lib/
	else
                fftw-inc = /home/hbabak/fftw-2.1.5_gcc/include/
                fftw-lib = /home/hbabak/fftw-2.1.5_gcc/lib/
	endif
endif

ifneq "$(findstring brutus,$(shell hostname))" ""
	mpi-inc = ${MPI_ROOT}/include/
	mpi-lib = ${MPI_ROOT}/lib/
	# fftw-inc = /cluster/work/infk/hbabak/apps/fftw-2.1.5/include/
	# fftw-lib = /cluster/work/infk/hbabak/apps/fftw-2.1.5/lib/
endif

ifneq "$(findstring rosa,$(shell hostname))" ""
	mpi-inc = $(CRAY_MPICH2_DIR)/include/
	mpi-lib = $(CRAY_MPICH2_DIR)/lib/
	fftw-inc = /users/petrosk/fftw-2.1.5_gcc/include/
	fftw-lib = /users/petrosk/fftw-2.1.5_gcc/lib/
endif

ifneq "$(findstring daint,$(shell hostname))" ""
        mpi-inc = $(CRAY_MPICH2_DIR)/include/
        mpi-lib = $(CRAY_MPICH2_DIR)/lib/
        fftw-inc = /users/petrosk/fftw-2.1.5_gcc/include/
        fftw-lib = /users/petrosk/fftw-2.1.5_gcc/lib/
endif

ifneq "$(findstring mira,$(shell hostname))" ""
        fftw-inc = /soft/libraries/alcf/current/xl/FFTW3/include/
        fftw-lib = /soft/libraries/alcf/current/xl/FFTW3/lib/
endif

ifneq "$(findstring cetus,$(shell hostname))" ""
        fftw-inc = /soft/libraries/alcf/current/xl/FFTW3/include/
        fftw-lib = /soft/libraries/alcf/current/xl/FFTW3/lib/
endif

ifneq "$(findstring fen,$(shell hostname))" ""
       fftw-inc = ${FFTW3_INC}
       fftw-lib = ${FFTW3_LIB}
endif

ifneq "$(findstring juq,$(shell hostname))" ""
        fftw-inc = ${FFTW3_INCLUDE}
        fftw-lib = ${FFTW3_LIB}
endif

CPPFLAGS += -I../../../Cubism/source/ -I../source -I../../MPCFcore/source/ -I../../MPCFnode/source/
CPPFLAGS +=  -I$(mpi-inc)

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

ifeq "$(fftw)"  "1"
#FFTW always in double precision
	CPPFLAGS += -I$(fftw-inc) -D_USE_FFTW_
	LIBS += -L$(fftw-lib) -lfftw3_omp -lfftw3_mpi -lfftw3
endif

##################
VPATH := ../source/ ../../MPCFcore/source/ ../../MPCFnode/source/ ../../../Cubism/source/
.DEFAULT_GOAL := mpcf-cluster

OBJECTS = ../../MPCFcore/makefiles/Convection_CPP.o \
		  ../../MPCFcore/makefiles/Source_CPP.o \
		  ../../MPCFcore/makefiles/Update.o \
		  ../../MPCFcore/makefiles/Update_State_CPP_5eq.o \
		  ../../MPCFcore/makefiles/MaxSpeedOfSound.o \
		  ../../MPCFcore/makefiles/MaxSpeedOfSound_CPP_5eq.o \
		  ../../MPCFcore/makefiles/Convection_CPP_HLLC_5eq.o \
		  ../../MPCFcore/makefiles/Convection_CPP_HLLC_5eq_nonuniform.o \
		  ../../MPCFcore/makefiles/InterfaceSharpening_CPP_5eq.o \
		  ../../MPCFcore/makefiles/InterfaceSharpening_CPP_cell_compact_5eq.o \
		  ../../MPCFcore/makefiles/MaxInterfaceVel_CPP_5eq.o \
		  ../../MPCFcore/makefiles/DivTensor_CPP.o \
		  ../../MPCFcore/makefiles/Diffusion_CPP.o \
		  ../../MPCFcore/makefiles/SurfaceTension_CPP.o
OBJECTS += ../../MPCFnode/makefiles/FlowStep_LSRK3.o \
		   ../../MPCFnode/makefiles/WaveletCompressor.o \
		   ../../MPCFnode/makefiles/Types.o \
		   ../../MPCFnode/makefiles/VectorOperator.o \
		   ../../MPCFnode/makefiles/Streamer.o \
		   ../../MPCFnode/makefiles/TestLabs.o \
		   ../../MPCFnode/makefiles/BoundaryGhostCells.o

OBJECTS += Profiler.o Histogram.o
OBJECTS += main.o \
		   FlowStep_LSRK3MPI.o

ifeq "$(fftw)" "1"
	OBJECTS += PoissonSolverScalarFFTW_MPI.o
endif

ifeq "$(qpx)" "1"
	OBJECTS += ../../MPCFcore/makefiles/WenoSOA2D_QPX.o
	OBJECTS += ../../MPCFcore/makefiles/DivSOA2D_QPX.o
	OBJECTS += ../../MPCFcore/makefiles/HLLCSOA2D_QPX_5eq.o
endif

ifeq "$(qpxemu)" "1"
	OBJECTS += ../../MPCFcore/makefiles/WenoSOA2D_QPX.o
	OBJECTS += ../../MPCFcore/makefiles/DivSOA2D_QPX.o
	OBJECTS += ../../MPCFcore/makefiles/HLLCSOA2D_QPX_5eq.o
endif

all: mpcf-cluster

mpcf-cluster: $(OBJECTS)
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

../../MPCFcore/makefiles/%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) $(extra) -c $^ -o $@

../../MPCFnode/makefiles/%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) $(extra) -c $^ -o $@

../../MPCFnode/makefiles/WaveletCompressor.o: WaveletCompressor.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) $(extra)  -c $^ -o $@

clean:
	rm -f *.o mpcf-cluster

cleanall:
	(cd ../../MPCFnode/makefiles; make cleanall)
	make clean;

