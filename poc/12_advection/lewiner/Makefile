.POSIX:
LINK = $(CC)
CFLAGS = -g -O2

M = main
O = \
main.o\
mc.o\

$M: $O; $(LINK) $O $(LDFLAGS) -o $@
.c.o:; $(CC) $(CFLAGS) $< -c -o $@
lint:; make CFLAGS='-g -Wall -Wextra -pedantic'

test: $M; atest test/*

.PHONY: clean lint dep test
clean:; rm -f -- $M $O

#dep:; ./tool/dep
#include dep.make
