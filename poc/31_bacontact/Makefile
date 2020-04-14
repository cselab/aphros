.POSIX:
.PHONY: clean run all u

CC99 = mpicc
CFLAGS = -O3 -g
Q_FLAGS = -I$(CHPREFIX)/include -lm

M = \
sessile3D\

all: $M
.c:
	. ba.setenv && \
	CC99=$(CC99) qcc $< $(CFLAGS) $(LDFLAGS) $(Q_FLAGS) -o $@

$M:  u
u:; mkdir -p .u && cp -r $(CHPREFIX)/include/io .u/io
run: all
	./sessile3D
clean:; rm -rf $M .u
