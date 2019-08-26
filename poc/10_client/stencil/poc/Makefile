.POSIX:
LINK = $(CC)
CFLAGS = -g -O2

M = main
O = main.o

main: $O; $(LINK) $O $(LDFLAGS) -o $@
.c.o:; $(CC) $(CFLAGS) $< -c -o $@
lint:; $(MAKE) CFLAGS='-Wall -Wextra -g'

.PHONY: clean lint
clean:; rm -f -- $M $O
