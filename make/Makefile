.POSIX:
.SUFFIXES:
.SUFFIXES: .o .cpp .c

L = libaphros_static.a
all: $L
include ../make/obj.mk

MPICXX = mpicxx
MPICC = mpicc
CXXFLAGS = -O2 -g
CFLAGS = -O2 -g
HDF_CFLAGS = `pkg-config --cflags hdf5-openmpi`
1 = -D_USE_AVX_=0 -D_USE_FPZIP_=0 -D_USE_AMGX_=0 -D_USE_HYPRE_=0 -D_USE_DIM1_=0 -D_USE_DIM2_=1 -D_USE_DIM3_=1 -D_USE_DIM4_=1

$L: $O; ar rv $@ $O && ranlib $@
.cpp.o:; $(MPICXX) $< -c $(CXXFLAGS) -I$(APHROS_PREFIX)/include $(HDF_CFLAGS) -I. $1 -o $@
.c.o:; $(MPICC) $< -c $(CFLAGS) -o $@

install: $L
	mkdir -p $(APHROS_PREFIX)/lib && cp $L $(APHROS_PREFIX)/lib/

include ../make/dep.mk
