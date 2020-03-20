include par.make

include $(shell ap.makesim)

par.conf par.make: gen par.py
	./gen

conf: par.conf

claen::
	rm -f par.make

cleandat::
	rm -vrf sep_????
