hdf = 0
vtk = 0
numa = 0
cvt = 0
##nthreads ?= 24

include ../../Makefile.config
#include ../Makefile.config.peh


#ifneq "$(CC)" "icc"
#	ifneq "$(config)" "release"
#		DIEGO_AVXGCC_ISSUE = -O1
#	endif
#endif
#CPPFLAGS += -D_ALIGNBYTES_=$(align) -D_BLOCKSIZE_=$(bs) -DWENOEPS=$(weps)
#CPPFLAGS += -I../../../Cubism/source

#CPPFLAGS += -D_L1PAR
CPPFLAGS += -D_L2PAR

###############

OBJECTS = main.o \
		  Convection_CPP.o \
		  Source_CPP.o \
		  Update.o \
		  Update_State_CPP_5eq.o \
		  MaxSpeedOfSound.o \
		  Diffusion_CPP.o \
		  DivTensor_CPP.o \
		  Convection_CPP_HLLC_5eq.o \
		  Convection_CPP_HLLC_5eq_nonuniform.o \
		  MaxSpeedOfSound_CPP_5eq.o \
		  MaxInterfaceVel_CPP_5eq.o \
		  InterfaceSharpening_CPP_5eq.o \
		  InterfaceSharpening_CPP_cell_compact_5eq.o
# OBJECTS = mymain.o Convection_CPP.o Source_CPP.o Update.o MaxSpeedOfSound.o Convection_CPP_HLLC.o Diffusion_CPP.o DivTensor_CPP.o Convection_CPP_HLLC_5eq.o MaxSpeedOfSound_CPP_5eq.o MaxInterfaceVel_CPP_5eq.o InterfaceSharpening_CPP_5eq.o InterfaceSharpening_CPP_cell_compact_5eq.o

ifeq "$(qpx)" "1"
OBJECTS += WenoSOA2D_QPX.o
OBJECTS += DivSOA2D_QPX.o
OBJECTS += HLLCSOA2D_QPX_5eq.o
NASTYFLAGS = -qlist -qsource -qreport -O5 -qarch=qp -qtune=qp  -qhot=level=2 -qunroll
endif

ifeq "$(qpxemu)" "1"
OBJECTS += WenoSOA2D_QPX.o
OBJECTS += DivSOA2D_QPX.o
OBJECTS += HLLCSOA2D_QPX_5eq.o
NASTYFLAGS = -Ofast
endif

VPATH := ../source/
.DEFAULT_GOAL := mpcf-core
header_files = $(wildcard ../source/*.h)

all: mpcf-core

mpcf-core: $(OBJECTS)
	$(CC) $(OPTFLAGS) $(CPPFLAGS) $^ -o $@ $(LIBS)

%.o: %.cpp $(header_files)
	$(CC) $(OPTFLAGS) $(CPPFLAGS) -c -o $@ $<

WenoSOA2D_QPX.o: WenoSOA2D_QPX.cpp WenoSOA2D_QPX.h common.h ../../MPCFthread/source/Weno_CPP.h ../../MPCFthread/source/Weno_QPX.h
	$(CC) $(NASTYFLAGS) $(CPPFLAGS) -c -o $@ $<

DivSOA2D_QPX.o: DivSOA2D_QPX.cpp DivSOA2D_QPX.h common.h
	$(CC) $(NASTYFLAGS) $(CPPFLAGS) -c -o $@ $<

HLLCSOA2D_QPX_5eq.o: HLLCSOA2D_QPX_5eq.cpp HLLCSOA2D_QPX_5eq.h common.h
	$(CC) $(NASTYFLAGS) $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f *.o mpcf-core
