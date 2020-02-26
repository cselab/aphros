#m = 32 32 1
m = 64 64 1
#m = 128 32 1
np = 1
#m = 512 128 1
#np = 256
#m = 256 64 1
#np = 64
bs = 32 32 1
tl = 1440

default:
	make cleanall
	make run

add.conf: add.conf.gen
	./$<

include $(shell ap.makesim)

conf: add.conf

clean::
	rm -vrf sc add.conf
