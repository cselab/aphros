m = 64 64 1
np = 4
bs = 32 32 1
tl = 1440

default:
	make cleanall
	make run

include $(shell ap.makesim)

conf: add.conf

clean::
	rm -vrf sc add.conf
