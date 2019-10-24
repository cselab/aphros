.POSIX:
.SUFFIXES: .o .c
CFLAGS = -g -Og
M_LDFLAGS = -lm
PREFIX = $(CHPREFIX)

O = main.o
L = libcsv.a
H = csv.h

libcsv.a: $O; ar rv $@ $O && ranlib $@
.o:; $(LINK) $< $(LDFLAGS) $(M_LDFLAGS) -o $@
.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
install: $L $H
	mkdir -p -- "$(PREFIX)/lib" "$(PREFIX)/include"
	cp -- $L "$(PREFIX)/lib"
	cp -- $H $(PREFIX)/inclde

.PHONY: clean lint install
clean:; rm -f -- $O $L
