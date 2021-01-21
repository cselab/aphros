.POSIX:
.SUFFIXES:
.SUFFIXES: .o .c
CFLAGS = -g -O2
GSL_CFLAGS = `gsl-config --cflags`
PREFIX = $(APHROS_PREFIX)
O = main.o err.o
L = libchmath.a
H = chmath.h
libchmath.a: $O; ar rv $@ $O && ranlib $@
main.o: chmath.h
main.o: err.h
err.o: err.h
.c.o:; $(CC) -c $(CFLAGS) $(GSL_CFLAGS) $<

lint:; make CFLAGS='-Wall -Wextra -g -O2'
install: $L $H
	mkdir -p -- "$(PREFIX)/lib" "$(PREFIX)/include/"
	cp -- $L "$(PREFIX)/lib/"
	cp -- $H "$(PREFIX)/include/"
clean:; rm -f -- $O $L
