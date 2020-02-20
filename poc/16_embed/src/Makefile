# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -fopenmp

all: check

madsen-sv.c: madsen.c
	ln -s madsen.c madsen-sv.c 
madsen-sv.tst: madsen.s
madsen-sv.tst: CFLAGS = -O2 -DSAINT_VENANT=1

madsen.tst: madsen-sv.tst
madsen.tst: CFLAGS = -O2

include $(BASILISK)/Makefile.defs
