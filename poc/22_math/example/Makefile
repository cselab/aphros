.POSIX:
M = vectors values
M_CFLAGS = -I$(CHPREFIX)/include `gsl-config --cflags`
M_LDFLAGS = -L$(CHPREFIX)/lib -lmath `gsl-config --libs`

LINK = $(CC)

all: $M
.o:
	$(LINK) $< $(LDFLAGS) $(M_LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $(M_CFLAGS) $< -c -o $@

.c:
%: %.c

test: $M; atest test/*
install:
lint:; make CFLAGS='-Wall -Wextra -g -O2'

.PHONY: clean test install run
clean:; rm -f $M $(M:=.o)
