include par.make

include $(shell ap.makesim)

par.conf par.make: gen par.py
	./gen | tee -a out

conf: par.conf add.conf

add.conf:
	touch $@

clean::
	rm -f par.make

cleandat::
	rm -vrf sep_????
