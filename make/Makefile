.POSIX:
.SUFFIXES:
.SUFFIXES: .o .cpp .c

L = libaphros.a
all: $L
include ../make/obj.mk

MPICXX = mpicxx
CXXFLAGS = -O2 -g
HDF_CFLAGS = `pkg-config --cflags hdf5-openmpi`
F = -D_USE_FPZIP_=0 -D_USE_AMGX_=0 -D_USE_HYPRE_=0

$L: $O; ar rv $@ $O && ranlib $@
.cpp.o:; $(MPICXX) $< -o $@ -c $(CXXFLAGS) -I$(APHROS_PREFIX)/include $(HDF_CFLAGS) -I. $F
.c.o:; $(MPICC) $< -c $(CFLAGS) -I. $F

include ../make/dep.mk
