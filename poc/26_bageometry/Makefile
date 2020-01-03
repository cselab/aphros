.POSIX:
.SUFFIXES: .o .c
CC = cc
CFLAGS = -g -Og
PREFIX = $(CHPREFIX)

O = main.o
L = libbageometry.a
H =

libbageometry.a: $O; ar rv $@ $O && ranlib $@

.c.o:; $(CC) -c $(CFLAGS) $<
main.o: bageometry.h def.h

.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
install: $L $H
	mkdir -p -- "$(PREFIX)/lib" "$(PREFIX)/include"
	cp -- $L "$(PREFIX)/lib"
	cp -- $H "$(PREFIX)/include"

.PHONY: clean lint install
clean:; rm -f -- $O $L
