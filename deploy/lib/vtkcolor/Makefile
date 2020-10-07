.POSIX:
CC = cc
LINK = $(CC)
CFLAGS = -g -Og
M_LDFLAGS = -lm

M = color

all: $M
.o:; $(LINK) $< $(LDFLAGS) $(M_LDFLAGS) -o $@
.c.o:; $(CC) $(CFLAGS) $< -c

.c:
%: %.c
.SUFFIXES: .o .c

install: all
	cp -v color ${APHROS_PREFIX}/bin/$(PRE)vtkcolor

lint:; make CFLAGS='-Wall -Wextra -g -O2'
.PHONY: clean lint
clean:; rm -f -- $M $(M:=.o)
