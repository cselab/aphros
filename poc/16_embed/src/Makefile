.SUFFIXES: .o .c
.PHONY: clean all
MDFLAGS = -fpic

include $(BASILISK)/config
QCC = $(BASILISK)/qcc

M = \
geometry.o\

all: $M
show:; echo $(CFLAGS)

.c.o:
	$(QCC) $(CFLAGS) -python -c $(MDFLAGS) -I$(PYTHONINCLUDE) $<

clean:; rm -f $M
