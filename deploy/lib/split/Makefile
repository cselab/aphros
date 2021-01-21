.POSIX:
PRE = ap.
LINK = $(CC)
CFLAGS = -g -Og
GSL_LDFLAGS = `gsl-config --libs`
M_LDFLAGS = -L$(APHROS_PREFIX)/lib -lvtk -lcsv -ltable -lchmath
M_CFLAGS = -I$(APHROS_PREFIX)/include

M = \
color\
column\
gyration\
filter\
split\

all: $M
.o:; $(LINK) $< $(LDFLAGS) $(M_LDFLAGS) $(GSL_LDFLAGS) -o $@
.c.o:; $(CC) $(CFLAGS) $(M_CFLAGS) $< -c

color.o: util.h
split.o: util.h

install: all
	@for i in $M; \
	do cp -v -- $$i $(APHROS_PREFIX)/bin/$(PRE)$$i; \
        done

.c:
%: %.c
%: %.sh
.SUFFIXES: .o .c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
.PHONY: clean lint install
clean:; rm -f -- $M $(M:=.o)
