# mesh size
m = 256 256 256
# number of processors
np = 256
# job time limit in minutes
tl = 600

include $(shell ap.makesim)

.PHONY: small

small:
	$(MAKE) cleanrun m='32 32 32' np=1
