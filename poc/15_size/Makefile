.POSIX:
LINK = $(CC)
CFLAGS = -g -Og
M_LDFLAGS = h.o -lm

M = \
area\
main\
rad\
size\
vtk2off\
vtk2off2\
hash\

all: $M
$M: h.o
.o:; $(LINK) $< $(LDFLAGS) $(M_LDFLAGS) -o $@
.c.o:; $(CC) $(CFLAGS) $< -c

.c:
%: %.c
.SUFFIXES: .o .c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
.PHONY: clean lint
clean:; rm -f -- $M $O
