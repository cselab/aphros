SHELL := /bin/bash

CC ?= mpic++
bs ?= 32
ap ?= float
config ?=release

include ../../Makefile.config

CPPFLAGS += -I../../Cubism/source -I../Cubism/source -I../../../Cubism/source/ -I../../MPCFcluster/source/ -I../../MPCFnode/source -I../../MPCcore/source -I.  -I../reader

##################
VPATH := ../../../Cubism/source/ ../../MPCFcluster/source/ ../../MPCFnode/source/ .
#.DEFAULT_GOAL := ch2ch ch2ch_io

OBJECTS = main.o WaveletCompressor.o
OBJECTS_IO = main_io.o WaveletCompressor.o

all: ch2ch ch2ch_io

ch2ch: $(OBJECTS) 
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

ch2ch_io: $(OBJECTS_IO) 
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

%.o: %.cpp 
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o ch2ch ch2ch_io
