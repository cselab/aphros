.POSIX:
LINK = $(CC)
CFLAGS = -g -O2
M_LDFLAGS = -lglut -lGL -lm

M = main
O = \
main.o\

$M: $O; $(LINK) $O $(LDFLAGS) $(M_LDFLAGS) -o $@
.c.o:; $(CC) $(CFLAGS) $< -c -o $@
lint:; make CFLAGS='-g -Wall -Wextra -pedantic'

.PHONY: clean lint
clean:; rm -f -- $M $O
