.POSIX:
.SUFFIXES: .o .c
CFLAGS = -g -Og
M_CFLAGS = `gsl-config --cflags`
PREFIX = $(CHPREFIX)

O = main.o err.o
L = libmath.a
H = math.h

libmath.a: $O; ar rv $@ $O && ranlib $@
main.o: math.h
main.o: err.h
err.o: err.h

.c.o:; $(CC) -c $(CFLAGS) $<
.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
install: $L $H
	mkdir -p -- "$(PREFIX)/lib" "$(PREFIX)/include/poc"
	cp -- $L "$(PREFIX)/lib/"
	cp -- $H "$(PREFIX)/include/poc/"

.PHONY: clean lint install
clean:; rm -f -- $O $L
