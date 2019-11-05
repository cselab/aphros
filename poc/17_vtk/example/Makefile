.POSIX:
LINK = $(CC)
CFLAGS = -g -Og
M_LDFLAGS = -L$(CHPREFIX)/lib -lvtk -lcsv -ltable -lm
M_CFLAGS = -I$(CHPREFIX)/include

M = \
add\
age\
field\
info\
ini\
off\
rad\
wall\
write\

all: $M
.o:; $(LINK) $< $(LDFLAGS) $(M_LDFLAGS) -o $@
.c.o:; $(CC) $(CFLAGS) $(M_CFLAGS) $< -c

age.o: util.h

.c:
%: %.c
%: %.sh
.SUFFIXES: .o .c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
.PHONY: clean lint
clean:; rm -f -- $M $O
