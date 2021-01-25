.POSIX:
.SUFFIXES:
.SUFFIXES: .c
CFLAGS = -g -O2
GSL_LDFLAGS = `gsl-config --libs`
CHMATH_FLAGS = -I.. -L.. -lchmath $(GSL_LDFLAGS)
M = eig_values
all: $M
.c:; $(CC) $< $(CFLAGS) $(CHMATH_FLAGS) $(LDFLAGS) -o $@
clean:; rm -f -- $M
