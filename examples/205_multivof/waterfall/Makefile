include par.make

include $(shell ap.makesim)

par.make par.conf: pargen par.py
	./pargen

par.py:
	touch "$@"

clean::
	rm -f _m _bs _np _tl par.conf b.dat par.make bc.dat
