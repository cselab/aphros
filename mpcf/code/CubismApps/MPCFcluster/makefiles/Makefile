SHELL := /bin/bash

CC ?= mpicxx
LD ?= mpicxx

include ../../Makefile.config

mpi-inc ?=.
mpi-lib ?=.

CPPFLAGS += -I../../../Cubism/source/ -I../source -I../../MPCFcore/source/ -I../../MPCFnode/source/
CPPFLAGS +=  -I$(mpi-inc)

ifeq "$(findstring mpi,$(CC))" ""
	LIBS += -L$(mpi-lib) -lmpi -lmpi_cxx
endif

##################
VPATH := ../source/ ../../MPCFcore/source/ ../../MPCFnode/source/ ../../../Cubism/source/
.DEFAULT_GOAL := mpcf-cluster

OBJECTS = main.o

all: mpcf-cluster

mpcf-cluster: $(OBJECTS)
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

../../MPCFcore/makefiles/%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) $(extra) -c $^ -o $@

../../MPCFnode/makefiles/%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) $(extra) -c $^ -o $@

clean:
	rm -f *.o mpcf-cluster

cleanall:
	(cd ../../MPCFnode/makefiles; make cleanall)
	make clean;

