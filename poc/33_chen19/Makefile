.POSIX:
.PHONY: clean all
.SUFFIXES: .off .scad

OPENSCAD = openscad

M = \
1.off\
2.off\

all: $M
.scad.off:
	$(OPENSCAD) $< -o $@

clean:; rm -f $M
