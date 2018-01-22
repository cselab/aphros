SHELL := /bin/bash

CC = mpic++
bs ?= 32
ap ?= float
hdf ?= 1
config ?=release

include ../../Makefile.config

nchannels ?= 1
CPPFLAGS += -DTOTAL_CHANNELS=$(nchannels)

# CSCS
#CPPFLAGS += -I${HDF5_INCLUDE_OPTS}

# BULLDOZER
#CPPFLAGS += -I${HOME}/usr/hdf5/include
#LIBS += -L${HOME}/usr/hdf5/lib -ldl

# PANDA
LIBS += -ldl

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


CPPFLAGS += -I../../Cubism/source -I../Cubism/source -I../../../Cubism/source/ -I../../MPCFcluster/source/ -I../../MPCFnode/source -I../../MPCcore/source -I.  -I../../tools/reader

##################
VPATH := ../../../Cubism/source/ ../../MPCFcluster/source/ ../../MPCFnode/source/ .
#.DEFAULT_GOAL := hdf2ch

OBJECTS = main.o WaveletCompressor.o

all: hdf2ch ch2diff

hdf2ch: $(OBJECTS) 
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

ch2diff: ch2diff.o WaveletCompressor.o
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

%.o: %.cpp 
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o hdf2ch ch2diff

cleanall:
	rm -f *.o hdf2ch out*
