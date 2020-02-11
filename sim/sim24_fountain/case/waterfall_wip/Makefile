include par.make

all:
	make cleanall
	make run

include $(shell ap.makesim)

par.make par.conf: pargen
	./pargen

clean::
	rm -f _m _bs _np _tl par.conf b.dat par.make
