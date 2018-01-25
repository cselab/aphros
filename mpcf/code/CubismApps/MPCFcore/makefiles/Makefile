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
		  Diffusion_CPP.o

ifeq "$(qpx)" "1"
NASTYFLAGS = -qlist -qsource -qreport -O5 -qarch=qp -qtune=qp  -qhot=level=2 -qunroll
endif

ifeq "$(qpxemu)" "1"
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

clean:
	rm -f *.o mpcf-core
