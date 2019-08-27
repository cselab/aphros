.POSIX:
CC = c99
FC = gfortran
LINK = $(CC)
FCFLAGS = -g -O2

M = main
O = \
main.o\
sz.o\

all: $M
$M: $O; $(LINK) $O $(LDFLAGS) -lm -o $@
.f90.o:; $(FC) $(FCFLAGS) $< -c -o $@
.c.o:; $(CC) $(CFLAGS) $< -c -o $@

.SUFFIXES: .o .f90 .c

lint:; $(MAKE) FCFLAGS='-Wall -Wextra -O2 -g'

.PHONY: clean lint
clean:; rm -f -- $O $M
