hdf = 0
vtk = 0
numa = 0
cvt = 0
##nthreads ?= 24

include ../../Makefile.config

ifeq "$(qpx)" "1"
NASTYFLAGS = -qsmp=omp:noauto -qlist -qsource -qreport -O5 -qarch=qp -qtune=qp  -qhot=level=2 -qunroll
#LIBS += /gpfs/DDNgpfs1/walkup/mpi_trace/bgq/lib/libmpihpm_smp.a \
#			/bgsys/drivers/ppcfloor/bgpm/lib/libbgpm.a \
#			-L/bgsys/drivers/ppcfloor/spi/lib -lSPI_upci_cnk
else
NASTYFLAGS = -O3 -fno-expensive-optimizations -falign-functions=16
endif

CPPFLAGS += -I../../MPCFcore/source -DNDEBUG


#ifneq "$(CC)" "icc"
#	ifneq "$(config)" "release"
#		DIEGO_AVXGCC_ISSUE = -O1
#	endif
#endif
#CPPFLAGS += -D_ALIGNBYTES_=$(align) -D_BLOCKSIZE_=$(bs) -DWENOEPS=$(weps)
#CPPFLAGS += -I../../../Cubism/source

###############

OBJECTS = main-weno.o 

VPATH := ../source/  ../../MPCFcore/source
.DEFAULT_GOAL := weno

all: weno

weno: $(OBJECTS)
	$(CC) $(CPPFLAGS) $(NASTYFLAGS) $(extra) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CC) $(CPPFLAGS) $(NASTYFLAGS) $(extra) -c $^ -o $@

clean:
	rm -f *.o weno
