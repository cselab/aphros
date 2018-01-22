CC = mpic++
bgq ?= 1
qpx = 1
bs ?= 32
#vs ?= 512
vs ?= 64
hpm = 0
config=release

include ../../Makefile.config

CPPFLAGS += -I../reader
CPPFLAGS += $(extra)
CPPFLAGS += $(OPTFLAGS)
# CSCS
#CPPFLAGS += -I${HDF5_INCLUDE_OPTS}

# BULLDOZER
#CPPFLAGS += -I${HOME}/usr/hdf5/include
#LIBS += -L${HOME}/usr/hdf5/lib -ldl

# JUQUEEN
#CPPFLAGS += -I/bgsys/local/hdf5/include
#LIBS += -L/bgsys/local/hdf5/lib -L/bgsys/local/szip/lib/ -lsz

# TUKEY
#CPPFLAGS += -I${HDF5_HOME}/include
#LIBS += -L${HDF5_HOME}/lib

#CPPFLAGS += -I${HOME}/usr_tukey/include
#LIBS += -L${HOME}/usr_tukey/lib

# CETUS
#CPPFLAGS += -I/soft/libraries/hdf5/current/cnk-xl/current/include
#LIBS += -L/soft/libraries/hdf5/current/cnk-xl/current/lib


#LIBS = $(LDFLAGS)
#ifeq "$(zlib)" "1"
#	MYFLAGS += -D_USE_ZLIB_
#endif
#
#ifeq "$(lz4)" "1"
#	MYFLAGS += -D_USE_LZ4_
#endif


#MYFLAGS =-D_USE_ZLIB_ -D_ALIGNBYTES_=16 -I../lz4 -mtune=native -march=native -ftree-vectorize -fstrict-aliasing  -D_BLOCKSIZE_=$(bs) -D_VOXELS_=$(vs) -D_FLOAT_PRECISION_ -D_SP_COMP_ -g -fopenmp -I../reader/ -I../../../Cubism/source/ -Ofast -D_QPXEMU_ 
CPPFLAGS += -D_VOXELS_=$(vs) -I../reader/ -I../../../Cubism/source/


vpcachier: main.o WaveletCompressor.o 
	$(CC) WaveletCompressor.o  main.o -o vpcachier  $(CPPFLAGS) $(LIBS) -lz -lstdc++

main.o: main.cpp
	$(CC) $(CPPFLAGS) -c main.cpp -o main.o
	
WaveletCompressor.o: ../../MPCFnode/source/WaveletCompressor.cpp
	$(CC) $(CPPFLAGS) -c ../../MPCFnode/source/WaveletCompressor.cpp -S
	$(CC) $(CPPFLAGS) -c ../../MPCFnode/source/WaveletCompressor.cpp -o WaveletCompressor.o
	
clean:
	rm vpcachier *.o
