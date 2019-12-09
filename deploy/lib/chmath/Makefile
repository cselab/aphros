.POSIX:
.SUFFIXES: .o .c
CFLAGS = -g -Og
GSL_CFLAGS = `gsl-config --cflags`
PREFIX = $(CHPREFIX)

O = main.o err.o
L = libchmath.a
H = chmath.h

libchmath.a: $O; ar rv $@ $O && ranlib $@
main.o: chmath.h
main.o: err.h
err.o: err.h

.c.o:; $(CC) -c $(CFLAGS) $(GSL_CFLAGS) $<
.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
install: $L $H
	mkdir -p -- "$(PREFIX)/lib" "$(PREFIX)/include/"
	cp -- $L "$(PREFIX)/lib/"
	cp -- $H "$(PREFIX)/include/"

.PHONY: clean lint install
clean:; rm -f -- $O $L
