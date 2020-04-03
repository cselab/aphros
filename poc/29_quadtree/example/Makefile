.POSIX:
.SUFFIXES: .o .c
include ../conf.mk
PREFIX = $(CHPREFIX)
M_CFLAGS = -I../lib
M_LDFLAGS = -L../lib -lquadtree

M = \
fig\

all: $M
.c.o:; $(CC) -c $(CFLAGS) $(M_CFLAGS) $<
.o:; $(LINK) $< $(LDFLAGS) $(M_LDFLAGS) -o $@
.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'

.PHONY: clean lint install
clean:; rm -f -- $M
