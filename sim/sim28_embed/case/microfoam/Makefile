include par.make

include $(shell ap.makesim)

par.conf par.make: gen par.py
	./gen

conf: par.conf add.conf

add.conf:
	touch $@

par.py:
	touch $@

clean::
	rm -f par.make par.conf b.dat bc.dat bg.dat body.dat
