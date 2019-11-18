.POSIX:
.SUFFIXES: .o .c
CFLAGS = -g -Og
M_LDFLAGS = -lm
PREFIX = $(CHPREFIX)

O = main.o err.o hash.o memory.o
L = libhe.a
H = he.h

libhe.a: $O; ar rv $@ $O && ranlib $@

err.o: err.h
hash.o: err.h
hash.o: hash.h
hash.o: memory.h
main.o: err.h
main.o: hash.h
main.o: he.h
main.o: memory.h

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
