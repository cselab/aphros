include par.make

include $(shell ch.makesim)

par.make par.conf bc.conf: pargen
	./pargen

clean::
	rm -f _m _bs _np _tl par.conf b.dat par.make
