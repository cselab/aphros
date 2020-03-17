.POSIX:
.SUFFIXES: .o .c
CFLAGS = -g -Og
PREFIX = $(CHPREFIX)

O = adamar.o
L = libadamar.a
H = adamar.h

libadamar.a: $O; ar rv $@ $O && ranlib $@
adamar.o: adamar.h

.c.o:; $(CC) -c $(CFLAGS) $<
.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
install: $L $H
	mkdir -p -- "$(PREFIX)/lib" "$(PREFIX)/include"
	cp -- $L "$(PREFIX)/lib"
	cp -- $H "$(PREFIX)/include"

main: main.c $L
	$(CC) $(CFLAGS) $(LDFLAGS) -L. main.c -ladamar -lm -o main

.PHONY: clean lint install
clean:; rm -f -- $O $L
