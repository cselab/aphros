include par.make

include $(shell ap-makesim)

par.make par.conf bc.conf: pargen bodygen
	./pargen
	./bodygen

clean::
	rm -f _m _bs _np _tl par.conf par.make body.dat
