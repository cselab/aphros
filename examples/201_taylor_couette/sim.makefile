m = 32 32 1
#m = 64 64 1
#m = 128 128 1
np = 1
bs = 32 32 1
tl = 1440

#default:; @:
#.PHONY: default

default:
	make cleanall
	make run

include $(shell ap.makesim)

clean::
	rm -vrf sc
