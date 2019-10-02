.POSIX:
LINK = $(CC)
CFLAGS = -g -Og
M_LDFLAGS = -lm

M = color main

all: $M
.o:; $(LINK) $< $(LDFLAGS) $(M_LDFLAGS) -o $@
.c.o:; $(CC) $(CFLAGS) $< -c

.c:
%: %.c
.SUFFIXES: .o .c

install: all
	cp -v color ${CHPREFIX}/bin/$(PRE)vtkcolor



lint:; make CFLAGS='-Wall -Wextra -g -O2'
.PHONY: clean lint
clean:; rm -f -- $M $O
