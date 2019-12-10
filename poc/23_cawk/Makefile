.POSIX:
LINK = $(CC)
CFLAGS = -g -Og
GSL_LDFLAGS = `gsl-config --libs`
M_LDFLAGS = -L$(CHPREFIX)/lib -lcsv -ltable -lmath  -lm
M_CFLAGS = -I$(CHPREFIX)/include

M = \
rdf\

all: $M
.o:; $(LINK) $< $(LDFLAGS) $(M_LDFLAGS) $(GSL_LDFLAGS) -o $@
.c.o:; $(CC) $(CFLAGS) $(M_CFLAGS) $< -c

.c:
%: %.c
.SUFFIXES: .o .c

csv2sph.o: ico.inc
csv2sph.o: util.h

lint:; make CFLAGS='-Wall -Wextra -g -O2'
.PHONY: clean lint
clean:; rm -f -- $M $(M:=.o)
