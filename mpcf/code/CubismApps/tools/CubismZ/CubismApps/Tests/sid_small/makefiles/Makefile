SHELL := /bin/bash

CC = mpicxx
LD = mpicxx
bs ?= 32
ap ?= float
config=release

include ../../../Makefile.config

CPPFLAGS += -I../../../Cubism/source/ -I../../../Compressor/source/ -I../source

##################
VPATH := ../../../Cubism/source/ ../../../Compressor/source/ ../source/
.DEFAULT_GOAL := tests

OBJECTS = main.o WaveletCompressor.o

all: tests

tests: $(OBJECTS) 
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o tests

cleanall:
	rm -f *.o tests out*
