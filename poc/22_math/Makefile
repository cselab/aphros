.POSIX:
.SUFFIXES: .o .c
CFLAGS = -g -Og
M_LDFLAGS = -lm
PREFIX = $(CHPREFIX)

O = main.o
L = libmath.a
H = math.h

libmath.a: $O; ar rv $@ $O && ranlib $@
main.o: math.h
main.o: err.h
err.o: err.h

.o:; $(LINK) $< $(LDFLAGS) $(M_LDFLAGS) -o $@
.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
install: $L $H
	mkdir -p -- "$(PREFIX)/lib" "$(PREFIX)/include"
	cp -- $L "$(PREFIX)/lib"
	cp -- $H "$(PREFIX)/include"

.PHONY: clean lint install
clean:; rm -f -- $O $L
