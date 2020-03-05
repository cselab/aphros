.POSIX:
.PHONY: clean lint all
.SUFFIXES: .o .c

M = \
line_alpha\
embed_interpolate\

M_CFLAGS = -I../lib
M_LDFLAGS = -L../lib -l2 `python-config --libs`

LINK = $(CC)

all: $M
.o:
	$(LINK) $< $(LDFLAGS) $(M_LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $(M_CFLAGS) $< -c -o $@

%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'
clean:; rm -f $M $(M:=.o)
