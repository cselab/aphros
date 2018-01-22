include ../../Makefile.config

CC = gcc

#########

VPATH := ../source/ ../../MPCFcore/source/ ../../../Cubism/source/
.DEFAULT_GOAL := mpcf-node

# NODE
OBJECTS = FlowStep_LSRK3.o \
		  Types.o \
		  VectorOperator.o \
		  Streamer.o \
		  TestLabs.o \
		  BoundaryGhostCells.o

# CORE
OBJECTS += ../../MPCFcore/makefiles/Convection_CPP.o \
		   ../../MPCFcore/makefiles/Source_CPP.o \
		   ../../MPCFcore/makefiles/Update.o \
		   ../../MPCFcore/makefiles/Update_State_CPP_5eq.o \
		   ../../MPCFcore/makefiles/MaxSpeedOfSound.o \
		   ../../MPCFcore/makefiles/Convection_CPP_HLLC_5eq.o \
		   ../../MPCFcore/makefiles/Convection_CPP_HLLC_5eq_nonuniform.o \
		   ../../MPCFcore/makefiles/DivTensor_CPP.o \
		   ../../MPCFcore/makefiles/Diffusion_CPP.o \
		   ../../MPCFcore/makefiles/MaxSpeedOfSound_CPP_5eq.o \
		   ../../MPCFcore/makefiles/InterfaceSharpening_CPP_cell_compact_5eq.o \
		   ../../MPCFcore/makefiles/InterfaceSharpening_CPP_5eq.o \
		   ../../MPCFcore/makefiles/MaxInterfaceVel_CPP_5eq.o
# TODO: (fabianw@mavt.ethz.ch; Wed 04 Jan 2017 02:40:08 PM CET) we do not have
# surface tension here.  needs to be added once FlowStep_LSRK3 has been
# refactored

# QPX/SSE
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

#on brutus flowstep performs better without native arch
ifeq "$(findstring brutus,$(shell hostname))" ""
	FSDIEGOEXTRA=$(extra)
endif

OBJECTS += Profiler.o

# MAIN PROGRAMS
MPCF = main.o
TESTKERNEL = kernelTestUtil.o

all: mpcf-node

mpcf-node: $(MPCF) $(OBJECTS)
	$(CC) $(OPTFLAGS) $(CPPFLAGS) $^ -o $@ $(LIBS)

test-lab: ../source/Test_Lab.cpp ../source/Types.h ../source/Types.cpp
	$(CC) $(OPTFLAGS) $(CPPFLAGS) ../source/Test_Lab.cpp ../source/Types.cpp -o TestLab $(LIBS)

test-kernel: $(TESTKERNEL) $(OBJECTS)
	$(CC) $(OPTFLAGS) $(CPPFLAGS) $^ -o $@ $(LIBS)

WaveletCompressor.o: WaveletCompressor.cpp
	$(CC) $(OPTFLAGS) $(CPPFLAGS)  -c $^ -o $@

%.o: %.cpp
	$(CC) $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@
	# $(CC) -std=c++0x $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

../../MPCFcore/makefiles/%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o mpcf-node

cleanall:
	(cd ../../MPCFcore/makefiles; make clean)
	make clean;
