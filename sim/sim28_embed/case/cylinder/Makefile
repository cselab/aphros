m = 128 64 1
bs = 32 32 1
np = 1
tl = 1440

include $(shell ap.makesim)

par.conf: par.py gen
	./gen

conf: par.conf

cleandat::
	rm -vrf sep_????
