.POSIX:
.PHONY: clean run all

CC99 = mpicc
CFLAGS = -Og -g
INSIDE_FLAGS = -I$(CHPREFIX)/include  -L$(CHPREFIX)/lib -linside
Q_FLAGS = -L$(BASILISK)/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11 -lm

M = \
distance\

all: $M
.c:
	CC99=$(CC99) qcc $< $(CFLAGS) $(INSIDE_FLAGS) $(LDFLAGS) $(Q_FLAGS) -o $@

run: all
	./distance
clean:; rm -rf $M .u
