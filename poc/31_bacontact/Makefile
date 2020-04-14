.POSIX:
.PHONY: clean run all

CC99 = mpicc
CFLAGS = -O3 -g
Q_FLAGS = `pkg-config --cflags hdf5` -I$(CHPREFIX)/include -lm

M = \
sessile3D\

all: $M
.c:
	. ba.setenv && \
	CC99=$(CC99) qcc $< $(CFLAGS) $(LDFLAGS) $(Q_FLAGS) -o $@

$M:  .u/bah5.h
.u/bah5.h:
	mkdir -p .u && \
	cp $(CHPREFIX)/include/bah5.h .u/bah5.h

run: all
	./sessile3D
clean:; rm -rf $M .u
