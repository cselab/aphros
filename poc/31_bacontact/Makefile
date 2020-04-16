.POSIX:
.PHONY: clean run all

CC99 = mpicc
CFLAGS = -O3 -g
Q_FLAGS = -I$(CHPREFIX)/include -lm

M = \
sessile3D\

all: $M
.c:
	CC99=$(CC99) qcc $< $(CFLAGS) $(LDFLAGS) $(Q_FLAGS) -o $@

$M:  .u
.u:; ba.u $(CHPREFIX)/cmake/aphros/io.conf
run: all
	./sessile3D
clean:; rm -rf $M .u
