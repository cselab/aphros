.POSIX:
.SUFFIXES: .o .c
CC = cc
CFLAGS = -g -Og
PREFIX = $(CHPREFIX)

O = \
bbox.o\
err.o\
main.o\
memory.o\
off.o\
predicate.o\

L = libinside.a
H = inside.h

$L: $O; ar rv $@ $O && ranlib $@

.c.o:; $(CC) -c $(CFLAGS) $<
include dep.mk

.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
install: $L $H
	mkdir -p -- "$(PREFIX)/lib" "$(PREFIX)/include"
	cp -- $L "$(PREFIX)/lib"
	cp -- $H "$(PREFIX)/include"

.PHONY: clean lint install
clean:; rm -f -- $O $L
