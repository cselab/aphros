.SUFFIXES: .o .c
.PHONY: clean all
MDFLAGS = -fpic
include $(BASILISK)/config
QCC = $(BASILISK)/qcc
CC = c99
M_CFLAGS = `python-config --includes`

L = \
lib2.a\
lib3.a\

2 = 2.o
3 = 3.o
all: $L
lib2.a: $2; ar r $@ $2 && ranlib $@
lib3.a: $3; ar r $@ $3 && ranlib $@
.c.o:
	CC99=$(CC) $(QCC) $(CFLAGS) -python -c $(MDFLAGS) $(M_CFLAGS) $<
clean:
	rm -f $L $2 $3 $(2:.o=.i) $(3:.o=.i)
