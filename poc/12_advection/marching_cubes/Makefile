.POSIX:
LINK = $(CC)
CFLAGS = -g -O2
M_LDFLAGS = -lglut -lGL -lm

M = main
O = \
main.o\
sample.o\


$M: $O; $(LINK) $O $(LDFLAGS) $(M_LDFLAGS) -o $@
.c.o:; $(CC) $(CFLAGS) $< -c -o $@
lint:; make CFLAGS='-g -Wall -Wextra -pedantic'

sample.o: sample.h
main.o: sample.h

.PHONY: clean lint
clean:; rm -f -- $M $O
