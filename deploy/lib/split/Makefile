.POSIX:
PRE = ch.
LINK = $(CC)
CFLAGS = -g -Og
M_LDFLAGS = -L$(CHPREFIX)/lib -lvtk -lcsv -ltable -lm
M_CFLAGS = -I$(CHPREFIX)/include

M = \
color\
split\

all: $M
.o:; $(LINK) $< $(LDFLAGS) $(M_LDFLAGS) -o $@
.c.o:; $(CC) $(CFLAGS) $(M_CFLAGS) $< -c

color.o: util.h
split.o: util.h

install: all
	@for i in $M; \
	do cp -v -- $$i $(CHPREFIX)/bin/$(PRE)$$i; \
        done

.c:
%: %.c
%: %.sh
.SUFFIXES: .o .c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
.PHONY: clean lint install
clean:; rm -f -- $M $O
