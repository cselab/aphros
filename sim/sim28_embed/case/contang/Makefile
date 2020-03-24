m = 64 64 1
#m = 128 128 1
np = 4
bs = 32 32 1
tl = 1440

#m = 32 32 1
#bs = 16 16 1

default:
	make cleanall
	make run

add.conf: add.conf.gen
	./$<

include $(shell ap.makesim)

conf: add.conf

clean::
	rm -vrf sc add.conf
