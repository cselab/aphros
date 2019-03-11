M = main
S = ../../../../../src
L = ../../lib

$M: $M.o

CH_CFLAGS = -I$L -I$S
CH_LDFLAGS = -L$L -L$S/build/solver -L$S/build/util -lsegment

$M: $M.o; $(CXX) $< $(CH_LDFLAGS) -o $@
.c.o:; $(CC) -c $(CFLAGS) $(CH_CFLAGS) $< -o $@

.SUFFIXES: .c .o

.PHONY: clean
clean:; rm -f $M $M.o
