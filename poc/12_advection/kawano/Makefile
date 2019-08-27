.POSIX:
FC = gfortran
LINK = $(CC)
FCFLAGS = -g -O2

#M = sz.o
O = sz.o

all: $O
#$M: $O; $(LINK) $O $(LDFLAGS) -o $@
.f90.o:; $(FC) $(FCFLAGS) $< -c -o $@

.SUFFIXES: .o .f90

lint:; $(MAKE) FCFLAGS='-Wall -Wextra -O2 -g'

.PHONY: clean lint
clean:; rm -f -- $O
