.POSIX:
LINK = $(CC)
CFLAGS = -g -Og
M_LDFLAGS = -L$(APHROS_PREFIX)/lib -lvtk -lcsv -ltable -lm
M_CFLAGS = -I$(APHROS_PREFIX)/include

M = \
add\
age\
color\
field\
info\
ini\
off\
one\
rad\
remove\
wall\
write\
xyz\

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
clean:; rm -f -- $M $(M:=.o)
