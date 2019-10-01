.POSIX:
LINK = $(CC)
CFLAGS = -g -Og
M_LDFLAGS =

M = main
O = main.o

$M: $O; $(LINK) $O $(LDFLAGS) $(M_LDFLAGS) -o $@
.c.o:; $(CC) $(CFLAGS) $< -c -o $@
lint:; make CFLAGS='-g -Wall -Wextra -pedantic'

.PHONY: clean lint dep
clean:; rm -f -- $M $O
