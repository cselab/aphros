.POSIX:
.PHONY: clean run all

CC99 = mpicc
Q_FLAGS = -lm

M = \
sessile3D\

all: $M
.c:
	. ba.setenv && \
	CC99=$(CC99) qcc $< $(CFLAGS) $(LDFLAGS) $(Q_FLAGS) -o $@

run: all
	./sessile3D
clean:; rm -f $M
