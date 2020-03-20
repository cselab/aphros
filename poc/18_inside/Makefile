.POSIX:
.SUFFIXES: .o .c
CC = cc
CFLAGS = -g -Og
PREFIX = $(CHPREFIX)

O = \
bbox.o\
main.o\
err.o\
memory.o\
predicate.o\

L = libtable.a
H = inside.h

libtable.a: $O; ar rv $@ $O && ranlib $@

.c.o:; $(CC) -c $(CFLAGS) $<

main.o: err.h memory.h
err.o: err.h
memory.o: memory.h

.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
install: $L $H
	mkdir -p -- "$(PREFIX)/lib" "$(PREFIX)/include"
	cp -- $L "$(PREFIX)/lib"
	cp -- $H "$(PREFIX)/include"

.PHONY: clean lint install
clean:; rm -f -- $O $L
