.POSIX:
.SUFFIXES: .o .c
CFLAGS = -g -Og
PREFIX = $(APHROS_PREFIX)

O = main.o err.o memory.o
L = libcsv.a
H = csv.h

libcsv.a: $O; ar rv $@ $O && ranlib $@

main.o: csv.h err.h memory.h
err.o: err.h
memory.o: memory.h

.c.o:; $(CC) -c $(CFLAGS) $<
.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
install: $L $H
	mkdir -p -- "$(PREFIX)/lib" "$(PREFIX)/include"
	cp -- $L "$(PREFIX)/lib"
	cp -- $H "$(PREFIX)/include"

.PHONY: clean lint install
clean:; rm -f -- $O $L
