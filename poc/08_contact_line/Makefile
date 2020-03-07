.POSIX:
.SUFFIXES: .c .o
.PHONY: all lint clean

M_LDFLAGS = -lm
LINK = $(CC)
CFLAGS = -Og -g
M = \
main\

all: $M
.o:; $(LINK) $< $(M_LDFLAGS) $(LDFLAGS) -o $@
.c.o:; $(CC) $(CFLAGS) $< -c

.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -Og'
clean:; rm -f $M $(M:=.o)
