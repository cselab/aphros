M = main
S = ../../../src

$M: $M.o

CH_CFLAGS = -I../lib -I$S
#CH_LDFLAGS = -L$S/build/solver -L$S/build/util -lsolver -lnormal -lsuspender

$M: $M.o; $(CC) $< $(CH_LDFLAGS) -o $@
.c.o:; $(CC) -c $(CFLAGS) $(CH_CFLAGS) $< -o $@

.SUFFIXES: .c .o

.PHONY: clean
clean:; rm -f $M $M.o
