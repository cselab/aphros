.POSIX:
.SUFFIXES: .o .c .2d.o .3d.o
CFLAGS = -g -Og
PREFIX = $(CHPREFIX)

2 = geometry.2d.o
3 = geometry.3d.o
L = libbasilisk2d.a libbasilisk3d.a
H = basilisk.h

all: $L
libbasilisk2d.a: $O; ar rv $@ $2 && ranlib $@
libbasilisk3d.a: $O; ar rv $@ $3 && ranlib $@

libbasilisk2d.a: $2
libbasilisk3d.a: $3

.c.2d.o:; $(CC) -Ddimension=2 -c $(CFLAGS) $< -o $@
.c.3d.o:; $(CC) -Ddimension=3 -c $(CFLAGS) $< -o $@

.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
install: $L $H
	mkdir -p -- "$(PREFIX)/lib" "$(PREFIX)/include"
	cp -- $L "$(PREFIX)/lib"
	cp -- $H "$(PREFIX)/include"

.PHONY: clean lint install
clean:; rm -f -- $2 $3 $L
