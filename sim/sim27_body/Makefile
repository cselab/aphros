include par.make

include $(shell ch.makesim)

par.make par.conf bc.conf: pargen bcgen
	./pargen
	./bcgen

clean::
	rm -f _m _bs _np _tl par.conf b.dat par.make bc.conf
