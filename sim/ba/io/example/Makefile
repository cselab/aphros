.POSIX:
.PHONY: clean all
.SUFFIXES: .c

CC99 = mpicc
CFLAGS = -O3 -g
Q_FLAGS = -I$(CHPREFIX)/include -lm

M = \
write\
facet\

all: $M
.c:
	. ba.setenv && \
	CC99=$(CC99) qcc $< $(CFLAGS) $(LDFLAGS) $(Q_FLAGS) -o $@

$M: .u
.u:; ba.env ba.u $(CHPREFIX)/cmake/aphros/io.conf
clean:; rm -rf $M .u
