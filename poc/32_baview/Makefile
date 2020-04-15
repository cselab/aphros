.POSIX:
.PHONY: clean run all

CC99 = mpicc
CFLAGS = -O3 -g
Q_FLAGS = -L$(BASILISK)/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11 -lm

M = \
distance\

all: $M
.c:
	. ba.setenv && \
	CC99=$(CC99) qcc $< $(CFLAGS) $(LDFLAGS) $(Q_FLAGS) -o $@

run: all
	./distance
clean:; rm -rf $M .u
