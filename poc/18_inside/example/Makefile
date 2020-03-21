.POSIX:
.SUFFIXES: .o .c
CFLAGS = -g -Og
PREFIX = $(CHPREFIX)
LINK = $(CC)
M_CFLAGS = -I../lib
M_LDFLAGS = -L../lib -linside

M = \
off\
point\
point2\

all: $M
.c.o:; $(CC) -c $(CFLAGS) $(M_CFLAGS) $<
.o:; $(LINK) $< $(LDFLAGS) $(M_LDFLAGS) -o $@
.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'

.PHONY: clean lint install
clean:; rm -f -- $M
