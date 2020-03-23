.POSIX:
.SUFFIXES: .o .c
include ../conf.mk
PREFIX = $(CHPREFIX)

O = \
kdtree.o\

L = libkdtree.a
H = kdtree.h

$L: $O; ar rv $@ $O && ranlib $@

.c.o:; $(CC) -c $(CFLAGS) $< -o $@
kdtree.o: kdtree.h

.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
install: $L $H
	mkdir -p -- "$(PREFIX)/lib" "$(PREFIX)/include"
	cp -- $L "$(PREFIX)/lib"
	cp -- $H "$(PREFIX)/include"

dep:; for i in *.c; do $(CPP) -MM $$i; done > dep.mk

.PHONY: clean lint install dep
clean:; rm -f -- $O $L
