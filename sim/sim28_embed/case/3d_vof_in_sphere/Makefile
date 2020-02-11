#m = 32 32 32
m = 64 64 64
#m = 128 32 1
np = 64
#m = 512 128 1
#np = 256
#m = 256 64 1
#np = 64
bs = 16 16 16
tl = 1440

default:
	killall ap.mfer || true
	make cleanall
	make run

include $(shell ap.makesim)

clean::
	rm -vrf sc
