.POSIX:
.PHONY: clean all
.SUFFIXES: .stl .scad

OPENSCAD = openscad

M = \
1.stl\
2.stl\
3.stl\

all: $M
.scad.stl:
	$(OPENSCAD) $< -o $@

clean:; rm -f $M
