.POSIX:
.PHONY: clean all
.SUFFIXES: .off .scad

OPENSCAD = openscad

M = \
1.off\
2.off\
3.off\

all: $M
.scad.off:
	trap 'rm -f /tmp/t.$$$$.off; exit 2' 1 2 3 4 15 && \
	$(OPENSCAD) $< -o /tmp/t.$$$$.off && \
	sed 's,^OFF ,OFF\n,g' /tmp/t.$$$$.off > $@ && \
	rm -f /tmp/t.$$$$.off

clean:; rm -f $M
