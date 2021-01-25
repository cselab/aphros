include par.make

default: help
.PHONY: default

include $(shell ap.makesim)

par.conf par.make: gen par.py
	./gen

conf: par.conf add.conf $(hook)

par.py:
	touch $@

clean::
	rm -f par.make par.conf b.dat bc.dat bg.dat body.dat
