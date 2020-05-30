m = 32 32 32
np = 1
tl = 600

MAKE=make

default: nx32

include $(shell ap.makesim)

.PHONY: nx32 nx128 default

nx32:
	$(MAKE) cleanrun m='32 32 32' np=1

nx128:
	$(MAKE) cleanrun m='128 128 128' np=128
